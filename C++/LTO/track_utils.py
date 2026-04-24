import json
import numpy as np
import casadi as ca
from scipy.interpolate import splprep, splev


class Track:
    def __init__(self, json_path=None, N_res=1000):
        self.N_res = N_res
        self.track_length = 0
        self.s_grid = None
        self.kappa = None
        self.n_left = None
        self.n_right = None
        self.X_spline = None
        self.Y_spline = None
        self.phi_spline = None

        if json_path:
            self.load_from_json(json_path)
        else:
            self._set_default_track()

    def get_center_line(self, s):
        # s can be a scalar or a numpy array
        u = np.array(s) / self.track_length
        x, y = splev(u, self.X_spline)
        return x, y

    def get_heading(self, s):
        u = np.array(s) / self.track_length
        dx, dy = splev(u, self.X_spline, der=1)
        return np.arctan2(dy, dx)

    def _set_default_track(self):
        self.track_length = 100
        self.s_grid = np.linspace(0, self.track_length, self.N_res + 1)
        self.kappa = ca.interpolant(
            "kappa", "linear", [self.s_grid], np.zeros(self.N_res + 1)
        )
        self.n_left = ca.interpolant(
            "n_left", "linear", [self.s_grid], np.ones(self.N_res + 1) * 5.0
        )
        self.n_right = ca.interpolant(
            "n_right", "linear", [self.s_grid], -np.ones(self.N_res + 1) * 5.0
        )
        self.banking = ca.interpolant(
            "banking", "linear", [self.s_grid], np.zeros(self.N_res + 1)
        )

    def load_from_json(self, json_path):
        print(f"Loading track from {json_path}")
        with open(json_path, "r") as f:
            data = json.load(f)

        x = np.array(data["X"])
        y = np.array(data["Y"])
        xo = np.array(data["X_o"])
        yo = np.array(data["Y_o"])
        xi = np.array(data["X_i"])
        yi = np.array(data["Y_i"])

        # Ensure track is closed for splprep (per=True)
        tck_c, _ = splprep([x, y], s=0, per=True)
        tck_o, _ = splprep([xo, yo], s=0, per=True)
        tck_i, _ = splprep([xi, yi], s=0, per=True)

        # Metric length calculation
        u_fine = np.linspace(0, 1, 10000)
        x_fine, y_fine = splev(u_fine, tck_c)
        dx = np.diff(x_fine)
        dy = np.diff(y_fine)
        self.track_length = np.sum(np.sqrt(dx**2 + dy**2))

        # Resample for CasADi interpolants
        u_grid = np.linspace(0, 1, self.N_res + 1)
        self.s_grid = u_grid * self.track_length

        # Evaluate center line for coordinates and curvature
        xc, yc = splev(u_grid, tck_c)
        dxc, dyc = splev(u_grid, tck_c, der=1)
        ddxc, ddyc = splev(u_grid, tck_c, der=2)

        # Curvature formula
        curvature = (dxc * ddyc - dyc * ddxc) / (dxc**2 + dyc**2) ** 1.5

        # Boundary distances (Simplified projection)
        n_left_vals = []
        n_right_vals = []
        u_search = np.linspace(0, 1, 1000)
        xo_s, yo_s = splev(u_search, tck_o)
        xi_s, yi_s = splev(u_search, tck_i)

        for i in range(len(u_grid)):
            p = np.array([xc[i], yc[i]])
            dist_o = np.min(np.sqrt((xo_s - p[0]) ** 2 + (yo_s - p[1]) ** 2))
            dist_i = np.min(np.sqrt((xi_s - p[0]) ** 2 + (yi_s - p[1]) ** 2))
            n_left_vals.append(dist_o)
            n_right_vals.append(-dist_i)

        self.kappa = ca.interpolant("kappa", "linear", [self.s_grid], curvature)
        self.n_left = ca.interpolant("n_left", "linear", [self.s_grid], n_left_vals)
        self.n_right = ca.interpolant("n_right", "linear", [self.s_grid], n_right_vals)

        # Store splines for coordinate reconstruction if needed
        self.X_spline = tck_c
        self.Y_spline = tck_c  # tck_c is a tuple ([u], [x,y], k)

        print(f"Track loaded: {round(self.track_length, 2)}m, {self.N_res} points.")
