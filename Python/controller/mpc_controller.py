import numpy as np
import scipy as sp
from scipy.linalg import block_diag
from qpsolvers import solve_qp

import sympy as sy
from sympy.physics import mechanics
from scipy.signal import cont2discrete

class SimModel(object):
    def __init__(self, param=None, NX=None, NU=None):
        assert param is not None
        assert NX is not None
        assert NU is not None

        self.param = param
        self.NX = NX
        self.NU = NU

        self.jacA, self.jacB = self.genLinModel()
        self.force = self.genDynamicEquation()

    def SYMPY_rh_eq(self):
        raise NotImplementedError('SYMPY_rh_eq is not implemented')

    def genJacobian(self):
        MAT = self.SYMPY_rh_eq()
        q = sy.symbols('q:{0}'.format(self.NX))
        u = sy.symbols('u:{0}'.format(self.NU))
        return MAT.jacobian(q), MAT.jacobian(u)

    def genDynamicEquation(self):
        q = sy.symbols("q:{0}".format(self.NX))
        u = sy.symbols("u:{0}".format(self.NU))
        return sy.lambdify([q,u], self.SYMPY_rh_eq(), [{'atan':np.arctan, 'atan2':np.arctan2}, "numpy"])

    def genLinModel(self):
        q = sy.symbols("q:{0}".format(self.NX))
        u = sy.symbols("u:{0}".format(self.NU))
        A, B = self.genJacobian()
        return (sy.lambdify([q,u], np.squeeze(A), [{'atan':np.arctan, 'atan2':np.arctan2}, "numpy"]),
                sy.lambdify([q,u], np.squeeze(B), [{'atan':np.arctan, 'atan2':np.arctan2}, "numpy"]))

    def genDModel(self, x, dq, u, dT=0.1):
        vector = np.hstack((x, dq))
        f = self.force(vector, u).T.flatten()

        A_c = np.array(self.jacA(vector, u))
        B_c = np.array(self.jacB(vector, u))
        g_c = f - A_c@vector - B_c@u

        B = np.hstack((B_c, g_c.reshape((-1,1))))

        A_d, B_d, _, _, _ = cont2discrete((A_c, B, 0, 0), dT)
        g_d = B_d[:,self.NU]
        B_d = B_d[:,0:self.NU]

        return A_d, B_d, g_d

    def PredictForwardEuler(self, x, dq, u, dt):
        vector = np.hstack((x, dq))
        d_vector = self.force(vector, u).T.flatten()
        vector = vector + dt * d_vector
        return vector[0:3], vector[3:6]

class NMPC():
    def __init__(self, dT=0.02, time_horizon = 20,
                 H = None, J = None, q = None,
                 RH = None, RJ = None, r = None,
                 H_N = None, J_N = None, q_N = None,
                 dmodel = None,
                 G = None, h = None,
                 normalization_x = None,
                 normalization_u = None,
                 x_ubounds=[], x_lbounds=[],
                 u_ubounds=[], u_lbounds=[]):

        assert H != None
        assert J != None
        assert H_N != None
        assert J_N != None
        assert dmodel != None

        self.dT = dT
        self.time_horizon = time_horizon
        self.model = dmodel

        self.x_l = np.asarray(x_lbounds,dtype=np.float64)
        self.x_u = np.asarray(x_ubounds,dtype=np.float64)
        self.u_l = np.asarray(u_lbounds,dtype=np.float64)
        self.u_u = np.asarray(u_ubounds,dtype=np.float64)

        self.NX = self.x_u.shape[0]
        self.NU = self.u_u.shape[0]

        class StaticStageCost():
            def __init__(self, weight):
                self.weight = np.asarray(weight)
            def __call__(self, x_guess, u_guess, x_ref):
                return self.weight
        class StaticValueFunction():
            def __init__(self, weight):
                self.weight = np.asarray(weight)
            def __call__(self, x_guess, x_ref):
                return self.weight

        self.H = H if callable(H) else StaticStageCost(H)
        self.J = J if callable(J) else StaticStageCost(J)
        self.H_N = H_N if callable(H_N) else StaticValueFunction(H_N)
        self.J_N = J_N if callable(J_N) else StaticValueFunction(J_N)
        self.R_H = RH if callable(RH) else StaticStageCost(RH)
        self.R_J = RJ if callable(RJ) else StaticStageCost(RJ)

        if q is None:
            q = np.zeros(self.NX)
        self.q = q
        if r is None:
            r = np.zeros(self.NU)
        self.r = r

        if G is None:
            self.G = None
            self.h = None
        else:
            assert h is not None
            self.G = G if callable(H) else StaticStageCost(H)
            self.h = h if callable(J) else StaticStageCost(J)

        if normalization_x is not None:
            self.Norm = np.diag(normalization_x*(time_horizon+1) + normalization_u*time_horizon)
            self.Norm_inv = np.linalg.inv(self.Norm)
        else:
            self.Norm = None
            self.Norm_inv = None

    def iterate_NMPC(self, x_guess, u_guess, x_ref, verbose=False, warmstart=False):
        T = self.time_horizon
        X_DIM = self.NX*(T+1)
        U_DIM = self.NU*(T)

        P_Q_blocks = []
        q_q_blocks = []
        P_R_blocks = []
        q_r_blocks = []

        for k in range(T+1):
            if k==T:
                P_Q_blocks.append(self.H_N(x_guess[:,k], x_ref[:,k]))
                q_q_blocks.append(self.J_N(x_guess[:,k], x_ref[:,k])+self.q)
            else:
                P_Q_blocks.append(self.H(x_guess[:,k], u_guess[:,k], x_ref[:,k]))
                q_q_blocks.append(self.J(x_guess[:,k], u_guess[:,k], x_ref[:,k])+self.q)

                P_R_blocks.append(self.R_H(x_guess[:,k], u_guess[:,k], x_ref[:,k]))
                q_r_blocks.append(self.R_J(x_guess[:,k], u_guess[:,k], x_ref[:,k])+self.r)

        P = block_diag(*P_Q_blocks,*P_R_blocks)
        q = np.hstack(q_q_blocks+q_r_blocks)
        P = 0.5*(P.T+P)

        Ad, Bd, gd = zip(*[self.model(q[:3], q[3:], u, self.dT)
                           for q, u in zip(x_guess.T, u_guess.T)])

        A = block_diag(*Ad)
        B = block_diag(*Bd)
        b = np.hstack((
                x_guess[:,0],
                - np.hstack(gd)
            ))
        A = np.block([
                [np.eye(self.NX), np.zeros((self.NX, X_DIM + U_DIM - self.NX))],
                [A, np.zeros((X_DIM - self.NX, self.NX)), B]
            ])
        A -= np.block([
                [np.zeros((self.NX, X_DIM+U_DIM))],
                [
                    np.zeros((X_DIM-self.NX, self.NX)),
                    np.eye(X_DIM - self.NX),
                    np.zeros_like(B)
                ]
            ])

        ### Track Constratint
        G = [
                self.G(x_g, x_r)
                for x_g,x_r
                in zip(x_guess.T, x_ref.T)
            ]
        h_block = [
                    np.asarray(self.G(x_g, x_r))@x_g - np.asarray(self.h(x_g, x_r) )
                    for x_g, x_r
                    in zip(x_guess.T, x_ref.T)
                ]
        G = np.hstack([
                block_diag(*G),
                np.zeros((T+1, U_DIM))
            ])
        h = np.hstack(h_block)

        x_l = np.tile(self.x_l, T+1)
        # Set trust region
        x_l[6::self.NX] = -0.2 + x_guess[6]
        u_l = np.tile(self.u_l, T)
        x_l = np.hstack((x_l, u_l))

        x_u = np.tile(self.x_u, T+1)
        # Set trust region
        x_u[6::self.NX] = 0.2 + x_guess[6]
        u_u = np.tile(self.u_u, T)
        x_u = np.hstack((x_u, u_u))

        #print([x.shape for x in [P, q, G, h, A, b, x_l, x_u]])

        try:
            if self.Norm is None:
                ret = solve_qp(P, q, G, h, A, b, x_l, x_u, solver='osqp')
            else:
                init_val = self.Norm_inv@np.hstack((x_guess.T.ravel(), u_guess.T.ravel()))
                #print("Equation Const", np.all(A@init_val==b))
                #print("InEquation Const", np.all(G@init_val<=h))
                #print("Lower bound Const", np.all(init_val>=x_l))
                #print("Upper bound Const", np.all(init_val<=x_u))
                #print("")

                ret = solve_qp(self.Norm@P@self.Norm,
                               q@self.Norm,
                               G@self.Norm, h,
                               self.Norm_inv[:X_DIM,:X_DIM]@A@self.Norm, self.Norm_inv[:X_DIM,:X_DIM]@b,
                               self.Norm_inv@x_l, self.Norm_inv@x_u,
                               initvals=init_val, solver='osqp')
                if ret[0] is not None:
                    ret = self.Norm@ret
        except Exception as e:
            print(e)
            return np.zeros_like(x_guess), np.zeros_like(u_guess), None

        #if ret.dtype != np.object:
        if ret is not None:
            ret_x = ret[:X_DIM].reshape((-1, self.NX)).T
            ret_u = ret[X_DIM:].reshape((-1, self.NU)).T
            return ret_x, ret_u, 0.
        else:
            return x_guess, u_guess, None

