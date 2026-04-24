import casadi as ca
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import json
import os
import argparse
from track_utils import Track


class VehicleParams:
    def __init__(self, param_file=None):
        # Default parameters
        self.m = 650.0
        self.Iz = 450.0
        self.lf = 2.0
        self.lr = 2.0
        self.g = 9.81

        # Aerodynamics
        self.Cl = 2.5
        self.rho = 1.225
        self.S = 1.0
        self.aero_split_front = 0.5

        # Track width
        self.tf = 0.8
        self.tr = 0.8
        self.Cm1 = 8000.0
        self.Cm2 = 172.0
        self.Cr0 = 180.0
        self.Cr2 = 0.7
        self.CBf = 4000.0
        self.CBr = 2000.0

        # Tire (Magic Formula)
        self.Bf = 10.0
        self.Cf = -1.38
        self.Df = 1.6
        self.Ef = 0.0
        self.Br = 10.0
        self.Cr = -1.38
        self.Dr = 1.6
        self.Er = 0.0

        # Tire constraints / scaling
        self.E_long = 1.2
        self.E_eps = 1.05

        # Constraints / Extra
        self.h = 0.5
        self.car_l = 4.0
        self.car_w = 1.8

        # Cost Weights
        self.q_time = 10.0
        self.q_r = 0.01
        self.q_alpha_r = 10.0
        self.q_n = 0.001
        self.q_mu = 0.01
        self.r_delta = 0.1
        self.r_d = 0.01
        self.r_B = 0.01
        self.r_diff_split = 10.0
        self.r_Ddelta = 10.0
        self.r_Dd = 0.01
        self.r_DB = 0.01

        if param_file and os.path.exists(param_file):
            with open(param_file, "r") as f:
                data = json.load(f)
                for k, v in data.items():
                    if hasattr(self, k):
                        setattr(self, k, v)


class BicycleModel:
    def __init__(self, p: VehicleParams):
        self.p = p

    def get_forces(self, x, u, kappa_s, F_Nf=None, F_Nr=None):
        vx = x[3]
        vy = x[4]
        omega = x[5]
        delta = u[0]
        D = u[1]
        B = u[2]

        alpha_f = ca.atan((vy + omega * self.p.lf) / (vx + 1e-6)) - delta
        alpha_r = ca.atan((vy - omega * self.p.lr) / (vx + 1e-6))

        if F_Nf is None or F_Nr is None:
            F_aero = 0.5 * self.p.rho * self.p.S * self.p.Cl * vx**2
            F_Nf = (
                self.p.aero_split_front * F_aero
                + (self.p.lr / (self.p.lf + self.p.lr)) * self.p.m * self.p.g
            )
            F_Nr = (1.0 - self.p.aero_split_front) * F_aero + (
                self.p.lf / (self.p.lf + self.p.lr)
            ) * self.p.m * self.p.g

        def pacejka(alpha, B, C, D, E):
            return D * ca.sin(
                C * ca.atan(B * alpha - E * (B * alpha - ca.atan(B * alpha)))
            )

        Fy_f = F_Nf * pacejka(alpha_f, self.p.Bf, self.p.Cf, self.p.Df, self.p.Ef)
        Fy_r = F_Nr * pacejka(alpha_r, self.p.Br, self.p.Cr, self.p.Dr, self.p.Er)

        F_M = self.p.Cm1 * D - self.p.Cm2 * D * vx
        F_Bf = -self.p.CBf * B
        F_Br = -self.p.CBr * B
        F_Fric_const = -self.p.Cr0
        F_Fric_aero = -0.5 * self.p.rho * self.p.S * self.p.Cr2 * vx**2

        Fx_f = F_Bf + F_Fric_const * 0.5
        Fx_r = F_M + F_Br + F_Fric_const * 0.5
        F_fric = F_Fric_aero

        return Fx_f, Fx_r, Fy_f, Fy_r, F_fric, alpha_f, alpha_r, F_M, F_Bf, F_Br

    def get_dynamics(self, x, u, kappa_s, F_Nf=None, F_Nr=None):
        n = x[1]
        mu = x[2]
        vx = x[3]
        vy = x[4]
        omega = x[5]
        delta = u[0]
        Fx_f, Fx_r, Fy_f, Fy_r, F_fric, af, ar, fm, fbf, fbr = self.get_forces(
            x, u, kappa_s, F_Nf, F_Nr
        )
        prog = (1 - kappa_s * n) / (vx * ca.cos(mu) - vy * ca.sin(mu) + 1e-6)
        dx = ca.vertcat(
            0,
            prog * (vy * ca.cos(mu) + vx * ca.sin(mu)),
            prog * omega - kappa_s,
            prog
            * (1.0 / self.p.m)
            * (
                Fx_r
                + F_fric
                - Fy_f * ca.sin(delta)
                + Fx_f * ca.cos(delta)
                + self.p.m * vy * omega
            ),
            prog
            * (1.0 / self.p.m)
            * (
                Fy_r
                + Fx_f * ca.sin(delta)
                + Fy_f * ca.cos(delta)
                - self.p.m * vx * omega
            ),
            prog
            * (1.0 / self.p.Iz)
            * (
                -Fy_r * self.p.lr
                + (Fx_f * ca.sin(delta) + Fy_f * ca.cos(delta)) * self.p.lf
            ),
        )
        return dx, prog, af, ar, Fy_f, Fy_r, Fx_f, Fx_r, fm, F_fric


class FourWheelModel:
    def __init__(self, p: VehicleParams):
        self.p = p

    def get_forces(self, x, u, kappa_s, F_N=None):
        # x: [t_inv, n, mu, vx, vy, omega]
        # u: [delta, D, B, diff_split]
        vx = x[3]
        vy = x[4]
        omega = x[5]
        delta = u[0]
        D = u[1]
        B = u[2]
        diff_split = u[3] if u.shape[0] > 3 else 0.0

        L = self.p.lf + self.p.lr
        TW = self.p.tf + self.p.tr

        delta_l = delta - TW / (2 * L) * delta**2
        delta_r = delta + TW / (2 * L) * delta**2

        # Slip angles
        af_l = (
            ca.atan((vy + self.p.lf * omega) / (vx - self.p.tf * omega + 1e-6))
            - delta_l
        )
        af_r = (
            ca.atan((vy + self.p.lf * omega) / (vx + self.p.tf * omega + 1e-6))
            - delta_r
        )
        ar_l = ca.atan((vy - self.p.lr * omega) / (vx - self.p.tr * omega + 1e-6))
        ar_r = ca.atan((vy - self.p.lr * omega) / (vx + self.p.tr * omega + 1e-6))

        # Static loads

        front_frac = self.p.lr / L
        rear_frac = self.p.lf / L
        left_frac = self.p.tr / TW
        right_frac = self.p.tf / TW

        F_Nlf_s = self.p.m * self.p.g * front_frac * left_frac
        F_Nrf_s = self.p.m * self.p.g * front_frac * right_frac
        F_Nlr_s = self.p.m * self.p.g * rear_frac * left_frac
        F_Nrr_s = self.p.m * self.p.g * rear_frac * right_frac

        # Aero loads
        F_aero = 0.5 * self.p.rho * self.p.S * self.p.Cl * vx**2
        aero_f_l = F_aero * self.p.aero_split_front * 0.5
        aero_f_r = F_aero * self.p.aero_split_front * 0.5
        aero_r_l = F_aero * (1 - self.p.aero_split_front) * 0.5
        aero_r_r = F_aero * (1 - self.p.aero_split_front) * 0.5

        if F_N is None:
            # Simple approximation if no opti variables provided
            F_Nlf = aero_f_l + F_Nlf_s
            F_Nrf = aero_f_r + F_Nrf_s
            F_Nlr = aero_r_l + F_Nlr_s
            F_Nrr = aero_r_r + F_Nrr_s
        else:
            F_Nlf, F_Nrf, F_Nlr, F_Nrr = F_N[0], F_N[1], F_N[2], F_N[3]

        def pacejka(alpha, B, C, D, E):
            return D * ca.sin(
                C * ca.atan(B * alpha - E * (B * alpha - ca.atan(B * alpha)))
            )

        # Normalized lateral forces (mu_y)
        mu_yf_l = pacejka(af_l, self.p.Bf, self.p.Cf, self.p.Df, self.p.Ef)
        mu_yf_r = pacejka(af_r, self.p.Bf, self.p.Cf, self.p.Df, self.p.Ef)
        mu_yr_l = pacejka(ar_l, self.p.Br, self.p.Cr, self.p.Dr, self.p.Er)
        mu_yr_r = pacejka(ar_r, self.p.Br, self.p.Cr, self.p.Dr, self.p.Er)

        # Longitudinal forces
        F_M = self.p.Cm1 * D - self.p.Cm2 * D * vx
        F_MWlr = 0.5 * (1 - diff_split) * F_M - 0.5 * self.p.CBr * B
        F_MWrr = 0.5 * (1 + diff_split) * F_M - 0.5 * self.p.CBr * B
        F_MWlf = -0.5 * self.p.CBf * B
        F_MWrf = -0.5 * self.p.CBf * B

        F_fric = -self.p.Cr0 - 0.5 * self.p.rho * self.p.S * self.p.Cr2 * vx**2

        return (
            (F_Nlf, F_Nrf, F_Nlr, F_Nrr),
            (mu_yf_l, mu_yf_r, mu_yr_l, mu_yr_r),
            (F_MWlf, F_MWrf, F_MWlr, F_MWrr),
            F_fric,
            (af_l, af_r, ar_l, ar_r),
        )

    def get_dynamics(self, x, u, kappa_s, F_N=None):
        n = x[1]
        mu = x[2]
        vx = x[3]
        vy = x[4]
        omega = x[5]
        delta = u[0]

        Loads, Mus, Fx, F_fric, Alphas = self.get_forces(x, u, kappa_s, F_N)
        F_Nlf, F_Nrf, F_Nlr, F_Nrr = Loads
        mu_yf_l, mu_yf_r, mu_yr_l, mu_yr_r = Mus
        Fx_lf, Fx_rf, Fx_lr, Fx_rr = Fx
        af_l, af_r, ar_l, ar_r = Alphas

        prog = (1 - kappa_s * n) / (vx * ca.cos(mu) - vy * ca.sin(mu) + 1e-6)

        # Longitudinal dynamics
        ax_long = (
            (Fx_lf + Fx_rf) * ca.cos(delta)
            + Fx_lr
            + Fx_rr
            + F_fric
            - (F_Nlf * mu_yf_l + F_Nrf * mu_yf_r) * ca.sin(delta)
            + self.p.m * vy * omega
        )

        # Lateral dynamics
        ay_lat = (
            (F_Nlr * mu_yr_l + F_Nrr * mu_yr_r)
            + (F_Nlf * mu_yf_l + F_Nrf * mu_yf_r) * ca.cos(delta)
            + (Fx_lf + Fx_rf) * ca.sin(delta)
            - self.p.m * vx * omega
        )

        # Yaw dynamics
        domega = (
            (F_Nlf * mu_yf_l + F_Nrf * mu_yf_r) * ca.cos(delta) * self.p.lf
            + (F_Nlf * mu_yf_l - F_Nrf * mu_yf_r) * ca.sin(delta) * self.p.tf
            + (Fx_lf + Fx_rf) * ca.sin(delta) * self.p.lf
            + (Fx_lf - Fx_rf) * ca.cos(delta) * self.p.tf
            - (F_Nlr * mu_yr_l + F_Nrr * mu_yr_r) * self.p.lr
            + (Fx_lr - Fx_rr) * self.p.tr
        )

        dx = ca.vertcat(
            0,
            prog * (vy * ca.cos(mu) + vx * ca.sin(mu)),
            prog * omega - kappa_s,
            prog * (1.0 / self.p.m) * ax_long,
            prog * (1.0 / self.p.m) * ay_lat,
            prog * (1.0 / self.p.Iz) * domega,
        )

        return dx, prog, Alphas, Mus, Fx, F_fric, Loads


class LapTimeOptimizer:
    def __init__(
        self, params: VehicleParams, track: Track, N_opt=500, model_type="bicycle"
    ):
        self.p = params
        self.track = track
        self.N_opt = N_opt
        self.model_type = model_type
        if model_type == "bicycle":
            self.model = BicycleModel(params)
            self.nu = 3
            self.nnn = 2
        else:
            self.model = FourWheelModel(params)
            self.nu = 4
            self.nnn = 4

    def solve(self):
        ds = self.track.track_length / (self.N_opt + 1)
        opti = ca.Opti()
        X = opti.variable(6, self.N_opt + 1)
        U = opti.variable(self.nu, self.N_opt + 1)
        dU = opti.variable(self.nu, self.N_opt)
        F_N_v = opti.variable(self.nnn, self.N_opt)

        # Guess
        vg = 20.0
        X_g = np.zeros((6, self.N_opt + 1))
        X_g[0, :] = 1.0 / vg
        X_g[3, :] = vg
        X_g[4, :] = 0.01 * np.random.randn(self.N_opt + 1)
        X_g[5, :] = 0.01 * np.random.randn(self.N_opt + 1)
        opti.set_initial(X, X_g)
        opti.set_initial(F_N_v, self.p.m * self.p.g / self.nnn)

        time_c = 0
        reg_c = 0
        for i in range(self.N_opt):
            s = i * ds
            ks = float(self.track.kappa(s))
            X_curr = X[:, i]
            U_curr = U[:, i]

            if self.model_type == "bicycle":
                (
                    k1_dx_f,
                    k1_p,
                    k1_af,
                    k1_ar,
                    k1_ff,
                    k1_fr,
                    k1_fxf,
                    k1_fxr,
                    k1_fm,
                    k1_fric,
                ) = self.model.get_dynamics(
                    X_curr, U_curr, ks, F_N_v[0, i], F_N_v[1, i]
                )

                L = self.p.lf + self.p.lr
                ax_s = (
                    k1_fxr
                    + k1_fric
                    - k1_ff * ca.sin(U[0, i])
                    + k1_fxf * ca.cos(U[0, i])
                )
                F_a = 0.5 * self.p.rho * self.p.S * self.p.Cl * X[3, i] ** 2
                opti.subject_to(
                    F_N_v[0, i]
                    == self.p.aero_split_front * F_a
                    + (self.p.lr / L) * self.p.m * self.p.g
                    - (0.5 * self.p.h / L) * ax_s
                )
                opti.subject_to(
                    F_N_v[1, i]
                    == (1 - self.p.aero_split_front) * F_a
                    + (self.p.lf / L) * self.p.m * self.p.g
                    + (0.5 * self.p.h / L) * ax_s
                )
                opti.subject_to(F_N_v[0, i] >= 50)
                opti.subject_to(F_N_v[1, i] >= 50)

                # Friction circle
                opti.subject_to(
                    ca.sqrt((self.p.E_long * k1_fm) ** 2 + k1_fr**2)
                    <= (F_N_v[1, i] * self.p.Dr * self.p.E_long)
                )
                opti.subject_to(
                    ca.sqrt((self.p.E_long * k1_fxf) ** 2 + k1_ff**2)
                    <= (F_N_v[0, i] * self.p.Df * self.p.E_long)
                )
                reg_c += self.p.q_alpha_r * k1_ar**2

            else:  # four_wheel
                k1_dx_f, k1_p, Alphas, Mus, Fx, F_fric, Loads = self.model.get_dynamics(
                    X_curr, U_curr, ks, F_N_v[:, i]
                )
                af_l, af_r, ar_l, ar_r = Alphas
                mu_yf_l, mu_yf_r, mu_yr_l, mu_yr_r = Mus
                Fx_lf, Fx_rf, Fx_lr, Fx_rr = Fx
                F_Nlf, F_Nrf, F_Nlr, F_Nrr = Loads

                L = self.p.lf + self.p.lr
                TW = self.p.tf + self.p.tr
                # Simplified load transfer to match Julia's logic structure
                ax_s = (
                    (Fx_lf + Fx_rf) * ca.cos(U[0, i])
                    + Fx_lr
                    + Fx_rr
                    + F_fric
                    - (F_Nlf * mu_yf_l + F_Nrf * mu_yf_r) * ca.sin(U[0, i])
                )
                ay_s = (
                    (F_Nlr * mu_yr_l + F_Nrr * mu_yr_r)
                    + (F_Nlf * mu_yf_l + F_Nrf * mu_yf_r) * ca.cos(U[0, i])
                    + (Fx_lf + Fx_rf) * ca.sin(U[0, i])
                )

                F_a = 0.5 * self.p.rho * self.p.S * self.p.Cl * X[3, i] ** 2
                F_Nf_s = self.p.m * self.p.g * (self.p.lr / L)
                F_Nr_s = self.p.m * self.p.g * (self.p.lf / L)

                opti.subject_to(
                    F_N_v[0, i]
                    == 0.5 * F_a * self.p.aero_split_front
                    + 0.5 * F_Nf_s
                    - (0.5 * self.p.h / L) * ax_s
                    - (0.5 * self.p.h / TW) * ay_s
                )  # LF
                opti.subject_to(
                    F_N_v[1, i]
                    == 0.5 * F_a * self.p.aero_split_front
                    + 0.5 * F_Nf_s
                    - (0.5 * self.p.h / L) * ax_s
                    + (0.5 * self.p.h / TW) * ay_s
                )  # RF
                opti.subject_to(
                    F_N_v[2, i]
                    == 0.5 * F_a * (1 - self.p.aero_split_front)
                    + 0.5 * F_Nr_s
                    + (0.5 * self.p.h / L) * ax_s
                    - (0.5 * self.p.h / TW) * ay_s
                )  # LR
                opti.subject_to(
                    F_N_v[3, i]
                    == 0.5 * F_a * (1 - self.p.aero_split_front)
                    + 0.5 * F_Nr_s
                    + (0.5 * self.p.h / L) * ax_s
                    + (0.5 * self.p.h / TW) * ay_s
                )  # RR

                for j in range(4):
                    opti.subject_to(F_N_v[j, i] >= 50)

                # Friction ellipses
                opti.subject_to(
                    ca.sqrt((self.p.E_long * Fx_lf) ** 2 + (F_N_v[0, i] * mu_yf_l) ** 2)
                    <= (self.p.E_eps * F_N_v[0, i] * self.p.Df * self.p.E_long)
                )
                opti.subject_to(
                    ca.sqrt((self.p.E_long * Fx_rf) ** 2 + (F_N_v[1, i] * mu_yf_r) ** 2)
                    <= (self.p.E_eps * F_N_v[1, i] * self.p.Df * self.p.E_long)
                )
                opti.subject_to(
                    ca.sqrt((self.p.E_long * Fx_lr) ** 2 + (F_N_v[2, i] * mu_yr_l) ** 2)
                    <= (self.p.E_eps * F_N_v[2, i] * self.p.Dr * self.p.E_long)
                )
                opti.subject_to(
                    ca.sqrt((self.p.E_long * Fx_rr) ** 2 + (F_N_v[3, i] * mu_yr_r) ** 2)
                    <= (self.p.E_eps * F_N_v[3, i] * self.p.Dr * self.p.E_long)
                )

                reg_c += self.p.q_alpha_r * (ar_l**2 + ar_r**2)
                reg_c += self.p.r_diff_split * U[3, i] ** 2

            # Integration (Forward Euler)
            k1_d = k1_dx_f[1:]
            X_next = X[1:, i] + ds * k1_d

            opti.subject_to(X[0, i] == k1_p)
            opti.subject_to(X[1:, i + 1] == X_next)

            time_c += ds * X[0, i] * self.p.q_time
            reg_c += (
                self.p.q_r * X[5, i] ** 2
                + self.p.q_n * X[1, i] ** 2
                + self.p.q_mu * X[2, i] ** 2
            )
            opti.subject_to(dU[:, i] == (U[:, i + 1] - U[:, i]) / (ds * X[0, i] + 1e-6))
            reg_c += 1.0 * U[1, i] * U[2, i]
            reg_c += (
                self.p.r_Ddelta * dU[0, i] ** 2
                + self.p.r_Dd * dU[1, i] ** 2
                + self.p.r_delta * U[0, i] ** 2
                + self.p.r_d * U[1, i] ** 2
                + self.p.r_B * U[2, i] ** 2
            )

            nl = float(self.track.n_left(s))
            nr = float(self.track.n_right(s))
            opti.subject_to(
                X[1, i]
                >= nr
                + 0.5
                * (self.p.car_l * ca.sin(X[2, i]) + self.p.car_w * ca.cos(X[2, i]))
            )
            opti.subject_to(
                X[1, i]
                <= nl
                - 0.5
                * (self.p.car_l * ca.sin(X[2, i]) + self.p.car_w * ca.cos(X[2, i]))
            )
            opti.subject_to(
                X[1, i]
                >= nr
                + 0.5
                * (-self.p.car_l * ca.sin(X[2, i]) + self.p.car_w * ca.cos(X[2, i]))
            )
            opti.subject_to(
                X[1, i]
                <= nl
                - 0.5
                * (-self.p.car_l * ca.sin(X[2, i]) + self.p.car_w * ca.cos(X[2, i]))
            )

        opti.subject_to(X[:, 0] == X[:, -1])
        opti.subject_to(U[:, 0] == U[:, -1])
        opti.subject_to(opti.bounded(5.0, X[3, :], 95.0))
        opti.subject_to(opti.bounded(-0.5, U[0, :], 0.5))
        opti.subject_to(opti.bounded(0.0, U[1, :], 1.0))
        opti.subject_to(opti.bounded(0.0, U[2, :], 1.0))
        if self.model_type == "four_wheel":
            opti.subject_to(opti.bounded(-0.25, U[3, :], 0.25))

        opti.minimize(time_c + reg_c)
        opti.solver(
            "ipopt", {"ipopt.max_iter": 500, "ipopt.tol": 1e-4, "ipopt.print_level": 5}
        )
        try:
            sol = opti.solve()
            return sol.value(X), sol.value(U), sol.value(F_N_v)
        except Exception:
            print("Solver failed, using debug values.")
            return opti.debug.value(X), opti.debug.value(U), opti.debug.value(F_N_v)


class LTOResult:
    def __init__(self, X, U, FN, lto: LapTimeOptimizer, track: Track):
        self.X = np.atleast_2d(X)
        self.U = np.atleast_2d(U)
        self.FN = np.atleast_2d(FN)
        self.lto = lto
        self.track = track
        self.params = lto.p
        self.model_type = lto.model_type
        self.N_opt = lto.N_opt

        # Post-processed data
        self.af_l = []
        self.ar_l = []
        self.us_l = []
        self.us2_l = []
        self.loads_l = []
        self.s_p = np.linspace(0, track.track_length, self.N_opt + 1)
        self.ds = track.track_length / self.N_opt

        self._process_results()

    def _process_results(self):
        if self.X.size == 0:
            return

        for i in range(self.N_opt + 1):
            s = i * self.ds
            ks = float(self.track.kappa(s))
            idx = i if i < self.N_opt else 0

            if self.model_type == "bicycle":
                _, _, af, ar, ff, fr, fxf, fxr, fm, fric = self.lto.model.get_dynamics(
                    self.X[:, i], self.U[:, i], ks, self.FN[0, idx], self.FN[1, idx]
                )
                self.af_l.append(np.degrees(float(af)))
                self.ar_l.append(np.degrees(float(ar)))
                self.us_l.append(
                    np.sqrt(float(fm) ** 2 + float(fr) ** 2)
                    / (float(self.FN[1, idx] * self.params.Dr) + 1e-6)
                )
                self.loads_l.append([float(self.FN[0, idx]), float(self.FN[1, idx])])
            else:
                _, _, Alphas, Mus, Fx, F_fric, Loads = self.lto.model.get_dynamics(
                    self.X[:, i], self.U[:, i], ks, self.FN[:, idx]
                )
                self.af_l.append(np.degrees(float(0.5 * (Alphas[0] + Alphas[1]))))
                self.ar_l.append(np.degrees(float(0.5 * (Alphas[2] + Alphas[3]))))
                usage_lr = np.sqrt(
                    (self.params.E_long * float(Fx[2])) ** 2
                    + (float(Loads[2]) * float(Mus[2])) ** 2
                ) / (self.params.E_eps * float(Loads[2]) * self.params.Dr + 1e-6)
                usage_rr = np.sqrt(
                    (self.params.E_long * float(Fx[3])) ** 2
                    + (float(Loads[3]) * float(Mus[3])) ** 2
                ) / (self.params.E_eps * float(Loads[3]) * self.params.Dr + 1e-6)
                self.us_l.append(float(usage_lr))
                self.us2_l.append(float(usage_rr))
                self.loads_l.append([float(load) for load in Loads])

        self.loads_l = np.array(self.loads_l)
        self.dt = self.ds * self.X[0, :-1]
        self.dot_delta = (self.U[0, 1:] - self.U[0, :-1]) / (self.dt + 1e-6)
        self.dot_D = (self.U[1, 1:] - self.U[1, :-1]) / (self.dt + 1e-6)
        self.dot_B = (self.U[2, 1:] - self.U[2, :-1]) / (self.dt + 1e-6)
        self.s_mid = 0.5 * (self.s_p[1:] + self.s_p[:-1])

        # Coordinate reconstruction
        self.tx = []
        self.ty = []
        for i in range(self.N_opt + 1):
            sv = i * self.ds
            xc, yc = self.track.get_center_line(sv)
            psi = self.track.get_heading(sv)
            self.tx.append(xc - self.X[1, i] * np.sin(psi))
            self.ty.append(yc + self.X[1, i] * np.cos(psi))


class LTOPlotter:
    def __init__(self, result: LTOResult):
        self.res = result

    def plot(self, html_path="lto_report.html"):
        num_rows = 11 if self.res.model_type == "four_wheel" else 10
        titles = [
            "Map",
            "vx [m/s]",
            "vy [m/s]",
            "omega [rad/s]",
            "n [m]",
            "mu [rad]",
            "t_inv [s/m]",
            "delta [rad]",
            "D [-]",
            "B [-]",
            "dot_delta [rad/s]",
            "dot_D [1/s]",
            "dot_B [1/s]",
            "Normal Loads [N]",
            "Slip Angles [deg]",
            "Usage [-]",
        ]
        if self.res.model_type == "four_wheel":
            titles.append("Diff Split [-]")

        fig = make_subplots(
            rows=num_rows,
            cols=2,
            specs=[[{"rowspan": 2}, {}], [None, {}]] + [[{}, {}]] * (num_rows - 2),
            subplot_titles=titles,
        )

        # Track plotting
        s_f = np.linspace(0, self.res.track.track_length, 500)
        cx, cy = self.res.track.get_center_line(s_f)
        fig.add_trace(
            go.Scatter(x=cx, y=cy, name="Center", line=dict(color="gray", dash="dash")),
            1,
            1,
        )

        lx, ly, rx, ry = [], [], [], []
        for sv in s_f:
            xc, yc = self.res.track.get_center_line(sv)
            psi = self.res.track.get_heading(sv)
            nl = float(self.res.track.n_left(sv))
            nr = float(self.res.track.n_right(sv))
            lx.append(xc - nl * np.sin(psi))
            ly.append(yc + nl * np.cos(psi))
            rx.append(xc - nr * np.sin(psi))
            ry.append(yc + nr * np.cos(psi))

        fig.add_trace(
            go.Scatter(x=lx, y=ly, name="Left", line=dict(color="black")), 1, 1
        )
        fig.add_trace(
            go.Scatter(x=rx, y=ry, name="Right", line=dict(color="black")), 1, 1
        )
        fig.add_trace(
            go.Scatter(
                x=self.res.tx,
                y=self.res.ty,
                name="Race Line",
                line=dict(color="red", width=3),
            ),
            1,
            1,
        )

        # Traces
        fig.add_trace(go.Scatter(x=self.res.s_p, y=self.res.X[3, :], name="vx"), 1, 2)
        fig.add_trace(go.Scatter(x=self.res.s_p, y=self.res.X[4, :], name="vy"), 2, 2)
        fig.add_trace(
            go.Scatter(x=self.res.s_p, y=self.res.X[5, :], name="omega"), 3, 1
        )
        fig.add_trace(go.Scatter(x=self.res.s_p, y=self.res.X[1, :], name="n"), 3, 2)
        fig.add_trace(go.Scatter(x=self.res.s_p, y=self.res.X[2, :], name="mu"), 4, 1)
        fig.add_trace(
            go.Scatter(x=self.res.s_p, y=self.res.X[0, :], name="t_inv"), 4, 2
        )

        # Inputs
        fig.add_trace(
            go.Scatter(x=self.res.s_p, y=self.res.U[0, :], name="delta"), 5, 1
        )
        fig.add_trace(go.Scatter(x=self.res.s_p, y=self.res.U[1, :], name="D"), 5, 2)
        fig.add_trace(go.Scatter(x=self.res.s_p, y=self.res.U[2, :], name="B"), 6, 1)

        # Input Rates
        fig.add_trace(
            go.Scatter(x=self.res.s_mid, y=self.res.dot_delta, name="dot_delta"), 6, 2
        )
        fig.add_trace(
            go.Scatter(x=self.res.s_mid, y=self.res.dot_D, name="dot_D"), 7, 1
        )
        fig.add_trace(
            go.Scatter(x=self.res.s_mid, y=self.res.dot_B, name="dot_B"), 7, 2
        )

        # Loads
        if self.res.model_type == "bicycle":
            fig.add_trace(
                go.Scatter(x=self.res.s_p, y=self.res.loads_l[:, 0], name="FNf"), 8, 1
            )
            fig.add_trace(
                go.Scatter(x=self.res.s_p, y=self.res.loads_l[:, 1], name="FNr"), 8, 1
            )
        else:
            fig.add_trace(
                go.Scatter(x=self.res.s_p, y=self.res.loads_l[:, 0], name="FNlf"), 8, 1
            )
            fig.add_trace(
                go.Scatter(x=self.res.s_p, y=self.res.loads_l[:, 1], name="FNrf"), 8, 1
            )
            fig.add_trace(
                go.Scatter(x=self.res.s_p, y=self.res.loads_l[:, 2], name="FNlr"), 8, 1
            )
            fig.add_trace(
                go.Scatter(x=self.res.s_p, y=self.res.loads_l[:, 3], name="FNrr"), 8, 1
            )

        # Slip Angles
        fig.add_trace(go.Scatter(x=self.res.s_p, y=self.res.af_l, name="af_avg"), 8, 2)
        fig.add_trace(go.Scatter(x=self.res.s_p, y=self.res.ar_l, name="ar_avg"), 8, 2)

        # Friction Usage
        fig.add_trace(
            go.Scatter(x=self.res.s_p, y=self.res.us_l, name="Usage (L)"), 9, 1
        )
        if self.res.model_type == "four_wheel":
            fig.add_trace(
                go.Scatter(x=self.res.s_p, y=self.res.us2_l, name="Usage (R)"), 9, 1
            )
            fig.add_trace(
                go.Scatter(x=self.res.s_p, y=self.res.U[3, :], name="diff_split"), 9, 2
            )

        fig.update_layout(height=2500, template="plotly_white", showlegend=True)
        fig.update_yaxes(scaleanchor="x", scaleratio=1, row=1, col=1)
        fig.write_html(html_path)
        print(f"Report saved to {html_path}")


class LTOExporter:
    def __init__(self, result: LTOResult):
        self.res = result

    def export_json(self, output_path="track_lto.json"):
        tx_raw = np.array(self.res.tx)
        ty_raw = np.array(self.res.ty)

        # Cumulative arc length along race line
        dx_race = np.diff(tx_raw)
        dy_race = np.diff(ty_raw)
        ds_race = np.sqrt(dx_race**2 + dy_race**2)
        s_race_raw = np.concatenate(([0], np.cumsum(ds_race)))
        L_race = s_race_raw[-1]

        # Equidistant grid along race line
        N_res_lto = 1000
        s_race_equi = np.linspace(0, L_race, N_res_lto)

        # Interpolate states to equidistant race line grid
        s_ref_at_equi = np.interp(s_race_equi, s_race_raw, self.res.s_p)
        tx_equi = np.interp(s_race_equi, s_race_raw, tx_raw)
        ty_equi = np.interp(s_race_equi, s_race_raw, ty_raw)
        vx_equi = np.interp(s_race_equi, s_race_raw, self.res.X[3, :])
        n_equi = np.interp(s_race_equi, s_race_raw, self.res.X[1, :])

        lx_equi, ly_equi, rx_equi, ry_equi = [], [], [], []
        nl_equi, nr_equi = [], []

        for i in range(N_res_lto):
            sv = s_ref_at_equi[i]
            xc, yc = self.res.track.get_center_line(sv)
            psi = self.res.track.get_heading(sv)
            nl = float(self.res.track.n_left(sv))
            nr = float(self.res.track.n_right(sv))

            lx_equi.append(xc - nl * np.sin(psi))
            ly_equi.append(yc + nl * np.cos(psi))
            rx_equi.append(xc - nr * np.sin(psi))
            ry_equi.append(yc + nr * np.cos(psi))
            nl_equi.append(nl - n_equi[i])
            nr_equi.append(nr - n_equi[i])

        track_lto = {
            "X": tx_equi.tolist(),
            "Y": ty_equi.tolist(),
            "X_o": lx_equi,
            "Y_o": ly_equi,
            "X_i": rx_equi,
            "Y_i": ry_equi,
            "s": s_race_equi.tolist(),
            "n_left": nl_equi,
            "n_right": nr_equi,
            "velocity": vx_equi.tolist(),
        }

        with open(output_path, "w") as f:
            json.dump(track_lto, f)
        print(f"LTO track saved to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Lap Time Optimizer")
    parser.add_argument(
        "--model",
        type=str,
        default="bicycle",
        choices=["bicycle", "four_wheel"],
        help="Vehicle model type",
    )
    parser.add_argument(
        "--N", type=int, default=500, help="Number of optimization intervals"
    )
    parser.add_argument(
        "--export", action="store_true", default=False, help="Export result to JSON"
    )
    parser.add_argument(
        "--plot", action="store_true", default=True, help="Generate HTML report"
    )
    parser.add_argument("--track", type=str, help="Path to track.json")
    parser.add_argument("--params", type=str, help="Path to model.json")
    parser.add_argument("--output", type=str, help="Path for output track_lto.json")
    parser.add_argument("--report", type=str, help="Path for output lto_report.html")
    args = parser.parse_args()

    WORK_DIR = os.path.dirname(os.path.abspath(__file__))

    # Default paths
    track_path = (
        args.track
        if args.track
        else os.path.join(os.path.dirname(WORK_DIR), "Params", "track.json")
    )
    params_path = (
        args.params
        if args.params
        else os.path.join(os.path.dirname(WORK_DIR), "Params", "model.json")
    )
    output_path = (
        args.output
        if args.output
        else os.path.join(os.path.dirname(WORK_DIR), "Params", "track_lto.json")
    )
    report_path = (
        args.report if args.report else os.path.join(WORK_DIR, "lto_report.html")
    )

    # Load track and parameters
    track = Track(track_path)
    params = VehicleParams(param_file=params_path)

    # Run optimization
    lto = LapTimeOptimizer(params, track, N_opt=args.N, model_type=args.model)
    X_o, U_o, FN_o = lto.solve()

    # Process results
    result = LTOResult(X_o, U_o, FN_o, lto, track)

    if result.X.size > 0:
        # Plotting
        if args.plot:
            plotter = LTOPlotter(result)
            plotter.plot(report_path)

        # Exporting
        if args.export:
            exporter = LTOExporter(result)
            exporter.export_json(output_path)
    else:
        print("Optimization failed to produce results.")
