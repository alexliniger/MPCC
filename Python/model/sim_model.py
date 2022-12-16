import numpy as np
from numba import jit, int32, float32, double, cfunc
from numba.experimental import jitclass

spec = [
    ('x', double[:]), ('dq', double[:]), ('u', double[:]),
    ('m', double), ('Iz', double),
    ('lf', double), ('lr', double),
    ('Bf', double), ('Cf', double), ('Df', double),
    ('Br', double), ('Cr', double), ('Dr', double),
    ('Cr0', double), ('Cr2', double),
    ('Cm1', double), ('Cm2', double),
    ('iterNum', int32), ('sim_dt', double), ('control_dt', double),
    ('car_shape', double[:,:]),
]

@jitclass(spec)
class VehicleSimModel(object):
    def __init__(self,  m=0.041, Iz=27.8E-6,
                        lf=0.029, lr=0.033,
                        Bf=2.579, Cf=1.2, Df=0.192,
                        Br=3.3852, Cr=1.2691, Dr=0.1737,
                        Cm1=0.287, Cm2=0.0545,
                        Cr0= 0.0518,
                        Cr2=0.00035, scale=1.0, control_dt = 10.0, sim_dt=1.0):

        self.x = np.asfortranarray(np.zeros(3, dtype=np.float64))
        self.dq = np.zeros(3, dtype=np.float64)
        self.u = np.zeros(2, dtype=np.float64)

        self.m = m
        self.Iz= Iz

        self.lf = lf
        self.lr = lr

        self.Bf = Bf
        self.Cf = Cf
        self.Df = Df
        self.Br = Br
        self.Cr = Cr
        self.Dr = Dr

        self.Cr0 = Cr0
        self.Cr2 = Cr2

        self.Cm1 = Cm1
        self.Cm2 = Cm2

        car_l = (lf + lr)/2 * scale
        car_w = car_l/2
        self.car_shape = np.asfortranarray(np.array([ [ car_l, car_w, 1.],
                                    [ car_l,-car_w, 1.],
                                    [-car_l,-car_w, 1.],
                                    [-car_l, car_w, 1.],
                                    [ car_l, car_w, 1.]], dtype=np.float64))
        self.sim_dt     = sim_dt
        self.control_dt = control_dt
        self.iterNum = int(self.control_dt/self.sim_dt)

    @property
    def shape(self):
        shape = np.dot(self.car_shape,
                       np.asfortranarray(np.array([
                           [ np.cos(self.x[2]),  np.sin(self.x[2]), 0.],
                           [-np.sin(self.x[2]),  np.cos(self.x[2]), 0.],
                           [ self.x[0]        ,  self.x[1]        , 1.]], dtype=np.float64)))
        return shape[:,:2]

    def ODE_rh_eq(self, x, dq, u):
        return self.ODE_rh_eq_x(x, dq, u), self.ODE_rh_eq_dq(x, dq, u)

    def ODE_rh_eq_dq(self, x, dq, u):
        ddq = np.zeros_like(dq)

        alpha_f =  np.arctan2(dq[1] + self.lf * dq[2], dq[0]) - u[1]
        alpha_r =  np.arctan2(dq[1] - self.lr * dq[2], dq[0])

        Ffy = - self.Df * np.sin(self.Cf*np.arctan(self.Bf*alpha_f))
        Fry = - self.Dr * np.sin(self.Cr*np.arctan(self.Br*alpha_r))
        Frx = self.Cm1*u[0]-self.Cm2*u[0]*dq[0]-self.Cr0-self.Cr2*dq[0]**2

        ddq[0] =  dq[1] * dq[2] + 1.0/self.m * (Frx - Ffy*np.sin(u[1]) )
        ddq[1] = -dq[0] * dq[2] + 1.0/self.m * (Fry + Ffy*np.cos(u[1]) )
        ddq[2] = 1.0/self.Iz * (Ffy * self.lf * np.cos(u[1]) - Fry * self.lr)

        return ddq

    def ODE_rh_eq_x(self, x, dq, u):
        dx = np.zeros_like(x)
        dx[0:2] = np.dot(np.array([[np.cos(x[2]), -np.sin(x[2])],
                                   [np.sin(x[2]),  np.cos(x[2])]]),
                         dq[0:2])
        dx[2] = dq[2]
        return dx

    def RK4(self, u, dt):
        k1x, k1dq = self.ODE_rh_eq(self.x, self.dq, u)
        k2x, k2dq = self.ODE_rh_eq(self.x+k1x*dt/2, self.dq+k1dq*dt/2, u)
        k3x, k3dq = self.ODE_rh_eq(self.x+k2x*dt/2, self.dq+k2dq*dt/2, u)
        k4x, k4dq = self.ODE_rh_eq(self.x+k3x*dt/2, self.dq+k3dq*dt/2, u)

        self.x  = self.x  + dt * (k1x/6 + k2x/3 + k3x/3 + k4x/6)
        self.dq = self.dq + dt * (k1dq/6 + k2dq/3 + k3dq/3 + k4dq/6)

    def LeapFrog(self, u, dt):
        dx = self.ODE_rh_eq_x(self.x, self.dq, u)
        ddq = self.ODE_rh_eq_dq(self.x, self.dq, u)
        self.x = self.x + dx * dt + ddq*dt**2/2

        ddq2= self.ODE_rh_eq_dq(self.x, self.dq, u)
        self.dq = self.dq + (ddq + ddq2)*dt/2

    def sim_step(self, u):
        for i in range(self.iterNum):
            self.RK4(u, self.sim_dt/1000)
            #self.LeapFrog(u, self.sim_dt/1000)

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    v = VehicleSimModel(scale=0.01)

    v.dq[0] = 0.2
    v.u[1] = np.deg2rad(20)

    x = np.array([v.x[0]])
    y = np.array([v.x[1]])
    car = []
    vhit = []

    for i in range(100):
        v.sim_step(v.u)
        x = np.append(x, v.x[0])
        y = np.append(y, v.x[1])

        carX, carY = v.shape.T
        car.append((carX, carY))
        vhit.append(v.dq[0])

    plt.figure(0)
    plt.plot(x, y, "ob-", label="trajectory")
    for carX, carY in car:
        plt.plot(carX, carY, 'k-')
    plt.axis("equal")
    plt.figure(1)
    plt.plot(range(len(vhit)), vhit)
    plt.grid(True)
    plt.show()
