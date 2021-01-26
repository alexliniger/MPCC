import numpy as np
from scipy import interpolate
from scipy import optimize

class CubicSplinePath:
    def __init__(self, x, y, close=False):
        x, y = map(np.asarray, (x, y))
        if close:
            x = np.append(x, x[0])
            y = np.append(y, y[0])
        s = np.append([0],np.cumsum((np.diff(x)**2 + np.diff(y)**2)**0.5))

        self.X = interpolate.CubicSpline(s, x, bc_type='periodic' if close else 'not-a-knot')
        self.Y = interpolate.CubicSpline(s, y, bc_type='periodic' if close else 'not-a-knot')

        self.dX = self.X.derivative(1)
        self.ddX = self.X.derivative(2)

        self.dY = self.Y.derivative(1)
        self.ddY = self.Y.derivative(2)

        self.length = s[-2] if close else s[-1]
    
    def calc_yaw(self, s):
        dx, dy = self.dX(s), self.dY(s)
        return np.arctan2(dy, dx)
    
    def calc_curvature(self, s):
        dx, dy   = self.dX(s), self.dY(s)
        ddx, ddy   = self.ddX(s), self.ddY(s)
        return (ddy * dx - ddx * dy) / ((dx ** 2 + dy ** 2)**(3 / 2))
    
    def __find_nearest_point(self, s0, x, y):
        def calc_distance(_s, *args):
            _x, _y= self.X(_s), self.Y(_s)
            return (_x - args[0])**2 + (_y - args[1])**2
        
        def calc_distance_jacobian(_s, *args):
            _x, _y = self.X(_s), self.Y(_s)
            _dx, _dy = self.dX(_s), self.dY(_s)
            return 2*_dx*(_x - args[0])+2*_dy*(_y-args[1])

        minimum = optimize.fmin_cg(calc_distance, s0, calc_distance_jacobian, args=(x, y), full_output=True, disp=False)
        return minimum

    def calc_track_error(self, x, y, s0):
        ret = self.__find_nearest_point(s0, x, y)
        
        s = ret[0][0]
        e = ret[1]

        k   = self.calc_curvature(s)
        yaw = self.calc_yaw(s)

        return e, k, yaw, s
