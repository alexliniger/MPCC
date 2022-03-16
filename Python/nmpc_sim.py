import json, os, argparse

import numpy as np
import scipy as sp
import pandas as pd

from util.cubic_spline_path import CubicSplinePath
from controller.mpc_controller import SimModel, NMPC
from model.sim_model import VehicleSimModel

base_path = "../C++/Params/"
def load_param(file_name):
    with open(os.path.join(base_path, file_name + ".json")) as json_file:
        return json.load(json_file)

class ExBicycleModel(SimModel):
    def SYMPY_rh_eq(self):
        import sympy as sy
        # Generalized coordinates and velocities
        ## Define q as [X,Y, phi, vx, vy, r, d, sigma, v_theta]
        q  = sy.symbols('q:{0}'.format(self.NX))
        dq = [q[3], q[4], q[5]]
        u  = [q[7], q[8], q[9]]
        ## Define du as [dd, dsigma, dv_theta]
        du  = sy.symbols('u:{0}'.format(self.NU))

        alpha_f =  sy.atan2(dq[1] + self.param['lf'] * dq[2], dq[0]) - u[1]
        alpha_r =  sy.atan2(dq[1] - self.param['lr'] * dq[2], dq[0])

        Ffy = - self.param['Df'] * sy.sin(self.param['Cf']*sy.atan(self.param['Bf']*alpha_f))
        Fry = - self.param['Dr'] * sy.sin(self.param['Cr']*sy.atan(self.param['Br']*alpha_r))
        Frx = self.param['Cm1'] * u[0] - self.param['Cm2']* u[0]*dq[0] - self.param['Cr0'] \
                 - self.param['Cr2'] * (dq[0] ** 2)

        F_MAT = sy.Matrix([
            sy.cos(q[2])*dq[0] - sy.sin(q[2])*dq[1],
            sy.sin(q[2])*dq[0] + sy.cos(q[2])*dq[1],
            dq[2],
            dq[1] * dq[2] + 1.0/self.param['m'] * (Frx - Ffy*sy.sin(u[1]) ),
           -dq[0] * dq[2] + 1.0/self.param['m'] * (Fry + Ffy*sy.cos(u[1]) ),
            1.0/self.param['Iz'] * (Ffy * self.param['lf'] * sy.cos(u[1]) - Fry * self.param['lr']),
            u[2],
            du[0], du[1], du[2]
        ])
        return F_MAT

    def PredictForwardEuler(self, x, u, dt):
        d_vector = self.force(x, u).T.flatten()
        return x + dt * d_vector

def main():
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--render",
        action="store_true",
        help="Render env states in a GUI window."
    )
    args = parser.parse_args()

    track = pd.read_json(os.path.join(base_path,"track.json"))
    track['S']=np.sqrt(track[['X', 'Y']].diff().pow(2).sum(axis=1)).cumsum()
    path_ref = CubicSplinePath(track['X'], track['Y'], close=True)

    dT = 0.02
    T  = 30
    simulation_time = 17

    bounds = load_param('bounds')
    cost   = load_param('cost')

    model_params = load_param('model')
    normalization = load_param('normalization')

    # Initialize model
    ## Define the environment (vehicle model)
    v = VehicleSimModel(scale=2, control_dt = dT*1000)
    ## Define controller internal model
    dmodel = ExBicycleModel(param=model_params, NX=10, NU=3)

    ### Theta state
    s0   = np.array([0.0, 1.0], dtype=np.float64)
    phi0 = np.arctan2(path_ref.dY(s0[0]), path_ref.dX(s0[0]))

    ### Vehicle state [X, Y, phi]
    x0   = np.array([path_ref.X(s0[0]), path_ref.Y(s0[0]), phi0],
                    dtype=np.float64)
    ### Vehicle state [vx, vy, r]
    dq0  = np.array([s0[1], 0., 0.],
                    dtype=np.float64)
    ### Vehicle input [d, sigma, v_theta]
    u0   = np.array([0., 0., s0[1]],
                    dtype=np.float64)

    ## Set the inital value to the internal variable of the simulation environment
    v.x = x0
    v.dq= dq0

    # Define Const Function
    def gen_cost_function():
        import sympy
        from sympy import cos, sin, atan2

        # Define x as [X,Y, phi, vx, vy, r, d, sigma, v_theta]
        x = sympy.symbols('x:{0}'.format(10))
        a = sympy.Function('X')(x[6])
        b = sympy.Function('Y')(x[6])

        xg = sympy.Matrix([_ for _ in x])
        xr = sympy.Matrix([a, b, 0, 0, 0, 0, 0, 0, 0, 0])

        X, Y, dX, dY, ddX, ddY = sympy.symbols('X Y dX dY ddX ddY')

        # [e_c; e_l] = U@(x-x_ref)
        U = sympy.Matrix([[ sin(atan2(b.diff(x[6]),a.diff(x[6]))),-cos(atan2(b.diff(x[6]),a.diff(x[6]))),0,0,0,0,0, 0,0,0],
                          [-cos(atan2(b.diff(x[6]),a.diff(x[6]))),-sin(atan2(b.diff(x[6]),a.diff(x[6]))),0,0,0,0,0, 0,0,0]])
        dE = sympy.derive_by_array(U*(xg-xr),x).reshape(10,2).transpose().tomatrix()
        dE = dE.subs(a.diff(x[6],2),ddX).subs(b.diff(x[6],2),ddY)
        dE = dE.subs(a.diff(x[6],1),dX).subs(b.diff(x[6],1),dY)
        dE = dE.subs(a,X).subs(b,Y)
        U  = U.subs(a.diff(x[6],1),dX).subs(b.diff(x[6],1),dY)

        U = sympy.lambdify([dX, dY], U, ["numpy"])
        dE = sympy.lambdify([x, X, Y, dX, dY, ddX, ddY], dE, ["numpy"])

        Q = np.diag([cost['qC'], cost['qL']])
        def H(x_guess, u_guess, x_ref):
            track_X = path_ref.X(x_guess[6])
            track_Y = path_ref.Y(x_guess[6])
            track_dX = path_ref.dX(x_guess[6])
            track_dY = path_ref.dY(x_guess[6])
            track_ddX = path_ref.ddX(x_guess[6])
            track_ddY = path_ref.ddY(x_guess[6])
            dE_val = dE(x_guess, track_X, track_Y, track_dX, track_dY, track_ddX, track_ddY)
            R = np.diag([0, 0, 0, 0, 0, 0, 0, cost['rD'], cost['rDelta'], cost['rVs']])
            return 2*dE_val.T@Q@dE_val + 2*R

        def J(x_guess, u_guess, x_ref):
            track_X = path_ref.X(x_guess[6])
            track_Y = path_ref.Y(x_guess[6])
            track_dX = path_ref.dX(x_guess[6])
            track_dY = path_ref.dY(x_guess[6])
            track_ddX = path_ref.ddX(x_guess[6])
            track_ddY = path_ref.ddY(x_guess[6])
            U_val = U(track_dX, track_dY)
            dE_val = dE(x_guess, track_X, track_Y, track_dX, track_dY, track_ddX, track_ddY)
            R = np.diag([0, 0, 0, 0, 0, 0, 0, cost['rD'], cost['rDelta'], cost['rVs']])
            return 2*((x_guess - x_ref).T@U_val.T@Q@dE_val -x_guess.T@dE_val.T@Q@dE_val).T + 2*x_guess@R

        def H_N(x_guess, x_ref):
            return 10*H(x_guess, None, x_ref)

        def J_N(x_guess, x_ref):
            return 10*J(x_guess, None, x_ref)

        q = np.zeros(10)
        q[-1] = -cost['qVs']

        R = 2*np.diag([cost['rdD'], cost['rdDelta'], cost['rdVs']])
        r = np.zeros(3)

        def R_J(x_guess, u_guess, x_ref):
            return 2*u_guess@R

        ## ieq constraint G+h<=0
        h = sympy.Matrix([(x[0] - a)**2 + (x[1] - b)**2 - model_params['R_in']**2])
        G = h.jacobian(x)

        G = G.subs(a.diff(x[6], 1), dX).subs(b.diff(x[6], 1), dY)
        G = G.subs(a,X).subs(b,Y)
        G = sympy.lambdify([x, X, Y, dX, dY], np.squeeze(G), ["numpy"])
        def G_func(x_guess, x_ref):
            track_X = path_ref.X(x_guess[6])
            track_Y = path_ref.Y(x_guess[6])
            track_dX = path_ref.dX(x_guess[6])
            track_dY = path_ref.dY(x_guess[6])
            return G(x_guess, track_X, track_Y, track_dX, track_dY)

        h = h.subs(a,X).subs(b,Y)
        h = sympy.lambdify([x, X, Y], np.squeeze(h), ["numpy"])
        def h_func(x_guess, x_ref):
            track_X = path_ref.X(x_guess[6])
            track_Y = path_ref.Y(x_guess[6])
            return h(x_guess, track_X, track_Y)

        return H, J, q, H_N, J_N, q, R, R_J, r, G_func, h_func

    H, J, q, H_N, J_N, _, RH, RJ, r, G, h = gen_cost_function()

    # Define NMPC Controller
    controller = NMPC(dT=dT, dmodel=dmodel.genDModel, time_horizon=T,
                  H=H, J=J, q=q, H_N=H_N, J_N=J_N,
                  RH= RH, RJ=RJ, r=r,
                  G=G, h=h,
                  #normalization_x = [normalization[e] for e in ["X", "Y", "phi", "vx", "vy", "r", "s", "D", "delta", "vs"]],
                  #normalization_u = [normalization[e] for e in ["dD", "dDelta", "dVs"]],
                  x_lbounds=[bounds['Xl'], bounds['Yl'],bounds['phil'],
                             bounds['vxl'],bounds['vyl'],bounds['rl'], bounds['sl'],
                             bounds['Dl'], bounds['deltal'], bounds['vsl']],
                  x_ubounds=[bounds['Xu'], bounds['Yu'],bounds['phiu'],
                             bounds['vxu'],bounds['vyu'],bounds['ru'], bounds['su'],
                             bounds['Du'], bounds['deltau'], bounds['vsu']],
                  u_lbounds=[bounds['dDl'], bounds['dDeltal'], bounds['dVsl']],
                  u_ubounds=[bounds['dDu'], bounds['dDeltau'], bounds['dVsu']])

    # Init Initial Guess
    u_bar = np.zeros((u0.shape[0], T), dtype=np.float64)
    x_bar = np.zeros((x0.shape[0]+1+dq0.shape[0]+u0.shape[0], T+1), dtype=np.float64)

    x_bar[0:3,0] = v.x
    x_bar[3:6,0] = v.dq
    x_bar[6,0]   = s0[0]
    x_bar[7:10,0]   = u0

    for t in range(T):
        x_bar[:,t+1] = dmodel.PredictForwardEuler(x_bar[:,t], u_bar[:,t], dT)

    # Logger
    x = []
    y = []
    car = []
    vlog = []

    render = args.render

    flag = False
    for i in range(int(simulation_time/dT)):
        e, _, _, s0[0] = path_ref.calc_track_error(v.x[0], v.x[1], s0[0])
        if e**0.5 > 0.15:
            break

        # Init or update the initial guess of SQP
        if flag:
            u_bar[:,0:-2] = u_bar[:,1:-1]
            x_bar[:,0:-2] = x_bar[:,1:-1]

            x_bar[0:3,0] = v.x
            x_bar[3:6,0] = v.dq
            x_bar[6, 0] = s0[0]
            x_bar[7:10, 0] = v.u

            x_bar[:,-2] = x_bar[:,-3]
            u_bar[:,-1] = np.zeros(3)
            x_bar[:,-1] = dmodel.PredictForwardEuler(x_bar[:,-2], u_bar[:,-1], dT)

            x_bar[6] = np.unwrap(x_bar[6], path_ref.length/2)
            x_bar[2] = np.unwrap(x_bar[2])

            updateW = 0.8
            QP_it = 3
        else:
            x_bar[0:3,0] = v.x
            x_bar[3:6,0] = v.dq
            x_bar[7:10,0] = u0
            u_bar = np.zeros((u0.shape[0], T), dtype=np.float64)

            x_bar[6] = np.unwrap(s0[0]+np.arange(0,T+1)*u0[2]*dT, path_ref.length/2)
            x_bar[2,1:] = np.arctan2(path_ref.dY(x_bar[6, 1:]), path_ref.dX(x_bar[6, 1:]))
            x_bar[2] = np.unwrap(x_bar[2])
            x_bar[0,1:] = path_ref.X(x_bar[6, 1:])
            x_bar[1,1:] = path_ref.Y(x_bar[6, 1:])
            x_bar[3,1:] = u0[2]

            updateW = 0.5
            QP_it = 5

            for t in range(T):
                x_bar[:,t+1] = dmodel.PredictForwardEuler(x_bar[:,t], u_bar[:,t], dT)

        # SQP Iteration
        flag = False
        for _ in range(QP_it):
            # MPC iteration
            x_ref = np.vstack([path_ref.X(x_bar[6]),  path_ref.Y(x_bar[6]),
                                np.arctan2(path_ref.dY(x_bar[6]), path_ref.dX(x_bar[6])),
                                np.zeros((3, T+1)),
                                x_bar[6],
                                np.zeros((3, T+1))])

            x_new, u_new, loss = controller.iterate_NMPC(x_bar, u_bar, x_ref)

            if loss is not None:
                x_bar  = updateW*x_bar + (1-updateW)*x_new
                u_bar  = updateW*u_bar + (1-updateW)*u_new
                flag = True

        # Get vehicle command 'u'
        v.u = x_bar[7:,1]

        x = np.append(x, v.x[0])
        y = np.append(y, v.x[1])

        carX, carY = v.shape.T
        car.append((carX, carY))
        vlog.append(np.hstack((v.x,v.dq)))

        if render:
            plt.cla()
            plt.plot(track['X'],  track['Y'], "r--" )
            plt.plot(track['X_i'],track['Y_i'], "k-")
            plt.plot(track['X_o'],track['Y_o'], "k-")

            plt.plot(carX,  carY, "b-")
            plt.plot(v.x[0], v.x[1], "x")
            plt.plot(x_ref[0], x_ref[1], "g--")
            plt.plot(x_bar[0], x_bar[1], "b-")
            plt.title("speed[m/s]:{:.2f}, deg:{:.2f}".format(v.dq[0], v.x[2]))
            plt.pause(0.001)

        # Step sim forward
        v.sim_step(v.u)

    # Plot result after the simulation
    plt.figure(0)
    plt.cla()
    plt.plot(track['X'],track['Y'], "r--")
    plt.plot(track['X_i'],track['Y_i'], "k-")
    plt.plot(track['X_o'],track['Y_o'], "k-")

    plt.plot(x, y, "b-", label="trajectory")
    #for carX, carY in car:
    #    plt.plot(carX, carY, 'k-')
    plt.axis("equal")
    plt.figure(1)
    plt.cla()
    plt.plot(range(len(vlog)), np.asarray(vlog)[:,3])
    plt.ylim(-0.5, 3.0)
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
