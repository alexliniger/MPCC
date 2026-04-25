// Copyright 2019 Alexander Liniger
// Licensed under the Apache License, Version 2.0

#include "Model/sim_model.h"

#include <algorithm>
#include <cmath>

namespace mpcc {

BicycleSimModel::BicycleSimModel(double Ts, const std::string& param_path)
    : Ts_(Ts), param_(param_path) {}

void BicycleSimModel::initialize(const State& x) {
  x_sim_.X = x.X;
  x_sim_.Y = x.Y;
  x_sim_.psi = x.phi;
  x_sim_.vx = x.vx;
  x_sim_.vy = x.vy;
  x_sim_.r = x.r;
  x_sim_.delta = x.delta;
  x_sim_.D = x.D;
  x_sim_.B = x.B;
}

void BicycleSimModel::step(const Input& input, double dt) {
  int n_sub = 10;
  double h = dt / n_sub;

  auto add_scaled_derivative = [](const SimState& s, const SimState& ds,
                                  double scale) {
    SimState res = s;
    res.X += ds.X * scale;
    res.Y += ds.Y * scale;
    res.psi += ds.psi * scale;
    res.vx += ds.vx * scale;
    res.vy += ds.vy * scale;
    res.r += ds.r * scale;
    res.delta += ds.delta * scale;
    res.D += ds.D * scale;
    res.B += ds.B * scale;
    return res;
  };

  for (int i = 0; i < n_sub; i++) {
    SimState k1 = getDerivative(x_sim_, input);
    SimState k2 =
        getDerivative(add_scaled_derivative(x_sim_, k1, h / 2.0), input);
    SimState k3 =
        getDerivative(add_scaled_derivative(x_sim_, k2, h / 2.0), input);
    SimState k4 = getDerivative(add_scaled_derivative(x_sim_, k3, h), input);

    x_sim_.X += (h / 6.0) * (k1.X + 2.0 * k2.X + 2.0 * k3.X + k4.X);
    x_sim_.Y += (h / 6.0) * (k1.Y + 2.0 * k2.Y + 2.0 * k3.Y + k4.Y);
    x_sim_.psi += (h / 6.0) * (k1.psi + 2.0 * k2.psi + 2.0 * k3.psi + k4.psi);
    x_sim_.vx += (h / 6.0) * (k1.vx + 2.0 * k2.vx + 2.0 * k3.vx + k4.vx);
    x_sim_.vy += (h / 6.0) * (k1.vy + 2.0 * k2.vy + 2.0 * k3.vy + k4.vy);
    x_sim_.r += (h / 6.0) * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r);
    x_sim_.delta +=
        (h / 6.0) * (k1.delta + 2.0 * k2.delta + 2.0 * k3.delta + k4.delta);
    x_sim_.D += (h / 6.0) * (k1.D + 2.0 * k2.D + 2.0 * k3.D + k4.D);
    x_sim_.B += (h / 6.0) * (k1.B + 2.0 * k2.B + 2.0 * k3.B + k4.B);
  }
}

BicycleSimModel::AeroForces BicycleSimModel::getAeroForces(
    const SimState& state) const {
  AeroForces aero;
  aero.F_downforce =
      0.5 * param_.rho * param_.S * param_.Cl * state.vx * state.vx;
  aero.F_drag = -0.5 * param_.rho * param_.S * param_.Cr2 * state.vx * state.vx;
  return aero;
}

BicycleSimModel::NormalForces BicycleSimModel::getNormalForces(
    const SimState& state, const AeroForces& aero) const {
  NormalForces F_N;
  F_N.F_N_front = param_.aero_split_front * aero.F_downforce +
                  (param_.lr / (param_.lf + param_.lr)) * param_.m * param_.g;
  F_N.F_N_rear = (1.0 - param_.aero_split_front) * aero.F_downforce +
                 (param_.lf / (param_.lf + param_.lr)) * param_.m * param_.g;
  return F_N;
}

BicycleSimModel::SlipAngles BicycleSimModel::getSlipAngles(
    const SimState& state) const {
  SlipAngles alpha;
  alpha.alpha_f =
      std::atan((state.vy + state.r * param_.lf) / std::max(state.vx, 1.0)) -
      state.delta;
  alpha.alpha_r =
      std::atan((state.vy - state.r * param_.lr) / std::max(state.vx, 1.0));
  return alpha;
}

BicycleSimModel::TireForces BicycleSimModel::getTireForces(
    const SimState& state, const SlipAngles& alpha,
    const NormalForces& F_N) const {
  TireForces F_tire;
  F_tire.Fy_f = F_N.F_N_front * param_.Df *
                std::sin(param_.Cf * std::atan(param_.Bf * alpha.alpha_f));
  F_tire.Fy_r = F_N.F_N_rear * param_.Dr *
                std::sin(param_.Cr * std::atan(param_.Br * alpha.alpha_r));

  F_tire.Fx_f = -param_.CBf * state.B - param_.Cr0 * 0.5;
  F_tire.Fx_r = param_.Cm1 * state.D - param_.Cm2 * state.D * state.vx -
                param_.CBr * state.B - param_.Cr0 * 0.5;
  return F_tire;
}

SimState BicycleSimModel::getDerivative(const SimState& s,
                                        const Input& u) const {
  SimState ds;

  AeroForces aero = getAeroForces(s);
  NormalForces F_N = getNormalForces(s, aero);
  SlipAngles alpha = getSlipAngles(s);
  TireForces F_tire = getTireForces(s, alpha, F_N);

  // Dynamics
  ds.X = s.vx * std::cos(s.psi) - s.vy * std::sin(s.psi);
  ds.Y = s.vx * std::sin(s.psi) + s.vy * std::cos(s.psi);
  ds.psi = s.r;
  ds.vx = 1.0 / param_.m *
          (F_tire.Fx_r + aero.F_drag - F_tire.Fy_f * std::sin(s.delta) +
           F_tire.Fx_f * std::cos(s.delta) + param_.m * s.vy * s.r);
  ds.vy = 1.0 / param_.m *
          (F_tire.Fy_r + F_tire.Fx_f * std::sin(s.delta) +
           F_tire.Fy_f * std::cos(s.delta) - param_.m * s.vx * s.r);
  ds.r = 1.0 / param_.Iz *
         (-F_tire.Fy_r * param_.lr +
          (F_tire.Fx_f * std::sin(s.delta) + F_tire.Fy_f * std::cos(s.delta)) *
              param_.lf);

  ds.delta = u.dDelta;
  ds.D = u.dD;
  ds.B = u.dB;

  return ds;
}

State BicycleSimModel::getMPCCState(double s_approx, double vs) const {
  State x;
  x.X = x_sim_.X;
  x.Y = x_sim_.Y;
  x.phi = x_sim_.psi;
  x.vx = x_sim_.vx;
  x.vy = x_sim_.vy;
  x.r = x_sim_.r;
  x.s = s_approx;
  x.D = x_sim_.D;
  x.B = x_sim_.B;
  x.delta = x_sim_.delta;
  x.vs = vs;
  return x;
}

// -----------------------------------------------------------------------------
// Four Wheel Model Implementation
// -----------------------------------------------------------------------------

FourWheelSimModel::FourWheelSimModel(double Ts, const std::string& param_path)
    : Ts_(Ts), param_(param_path) {}

void FourWheelSimModel::initialize(const State& x) {
  x_sim_.X = x.X;
  x_sim_.Y = x.Y;
  x_sim_.psi = x.phi;
  x_sim_.vx = x.vx;
  x_sim_.vy = x.vy;
  x_sim_.r = x.r;
  x_sim_.delta = x.delta;
  x_sim_.D = x.D;
  x_sim_.B = x.B;
}

void FourWheelSimModel::step(const Input& input, double dt) {
  int n_sub = 10;
  double h = dt / n_sub;

  auto add_scaled_derivative = [](const SimState& s, const SimState& ds,
                                  double scale) {
    SimState res = s;
    res.X += ds.X * scale;
    res.Y += ds.Y * scale;
    res.psi += ds.psi * scale;
    res.vx += ds.vx * scale;
    res.vy += ds.vy * scale;
    res.r += ds.r * scale;
    res.delta += ds.delta * scale;
    res.D += ds.D * scale;
    res.B += ds.B * scale;
    return res;
  };

  for (int i = 0; i < n_sub; i++) {
    SimState k1 = getDerivative(x_sim_, input);
    SimState k2 =
        getDerivative(add_scaled_derivative(x_sim_, k1, h / 2.0), input);
    SimState k3 =
        getDerivative(add_scaled_derivative(x_sim_, k2, h / 2.0), input);
    SimState k4 = getDerivative(add_scaled_derivative(x_sim_, k3, h), input);

    x_sim_.X += (h / 6.0) * (k1.X + 2.0 * k2.X + 2.0 * k3.X + k4.X);
    x_sim_.Y += (h / 6.0) * (k1.Y + 2.0 * k2.Y + 2.0 * k3.Y + k4.Y);
    x_sim_.psi += (h / 6.0) * (k1.psi + 2.0 * k2.psi + 2.0 * k3.psi + k4.psi);
    x_sim_.vx += (h / 6.0) * (k1.vx + 2.0 * k2.vx + 2.0 * k3.vx + k4.vx);
    x_sim_.vy += (h / 6.0) * (k1.vy + 2.0 * k2.vy + 2.0 * k3.vy + k4.vy);
    x_sim_.r += (h / 6.0) * (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r);
    x_sim_.delta +=
        (h / 6.0) * (k1.delta + 2.0 * k2.delta + 2.0 * k3.delta + k4.delta);
    x_sim_.D += (h / 6.0) * (k1.D + 2.0 * k2.D + 2.0 * k3.D + k4.D);
    x_sim_.B += (h / 6.0) * (k1.B + 2.0 * k2.B + 2.0 * k3.B + k4.B);
  }
}

FourWheelSimModel::AeroForces FourWheelSimModel::getAeroForces(
    const SimState& state) const {
  AeroForces aero;
  double F_aero = 0.5 * param_.rho * param_.S * param_.Cl * state.vx * state.vx;
  aero.F_downforce_f = F_aero * param_.aero_split_front * 0.5;
  aero.F_downforce_r = F_aero * (1.0 - param_.aero_split_front) * 0.5;
  aero.F_drag = -0.5 * param_.rho * param_.S * param_.Cr2 * state.vx * state.vx;
  return aero;
}

FourWheelSimModel::NormalForces FourWheelSimModel::getNormalForces(
    const SimState& state, const AeroForces& aero, double ax, double ay) const {
  NormalForces F_N;
  double L = param_.lf + param_.lr;
  double TW = param_.tf + param_.tr;

  double F_Nf_s = param_.m * param_.g * (param_.lr / L);
  double F_Nr_s = param_.m * param_.g * (param_.lf / L);

  F_N.F_N_fl = aero.F_downforce_f + 0.5 * F_Nf_s -
               (0.5 * param_.h_cg / L) * param_.m * ax -
               (0.5 * param_.h_cg / TW) * param_.m * ay;
  F_N.F_N_fr = aero.F_downforce_f + 0.5 * F_Nf_s -
               (0.5 * param_.h_cg / L) * param_.m * ax +
               (0.5 * param_.h_cg / TW) * param_.m * ay;
  F_N.F_N_rl = aero.F_downforce_r + 0.5 * F_Nr_s +
               (0.5 * param_.h_cg / L) * param_.m * ax -
               (0.5 * param_.h_cg / TW) * param_.m * ay;
  F_N.F_N_rr = aero.F_downforce_r + 0.5 * F_Nr_s +
               (0.5 * param_.h_cg / L) * param_.m * ax +
               (0.5 * param_.h_cg / TW) * param_.m * ay;

  // Prevent negative normal forces (lifting off)
  F_N.F_N_fl = std::max(F_N.F_N_fl, 0.0);
  F_N.F_N_fr = std::max(F_N.F_N_fr, 0.0);
  F_N.F_N_rl = std::max(F_N.F_N_rl, 0.0);
  F_N.F_N_rr = std::max(F_N.F_N_rr, 0.0);

  return F_N;
}

FourWheelSimModel::SlipAngles FourWheelSimModel::getSlipAngles(
    const SimState& state, double delta) const {
  SlipAngles alpha;
  double L = param_.lf + param_.lr;
  double TW = param_.tf + param_.tr;

  double delta_l = delta - (TW / (2.0 * L)) * delta * delta;
  double delta_r = delta + (TW / (2.0 * L)) * delta * delta;

  alpha.alpha_fl = std::atan((state.vy + param_.lf * state.r) /
                             std::max(state.vx - param_.tf * state.r, 1e-6)) -
                   delta_l;
  alpha.alpha_fr = std::atan((state.vy + param_.lf * state.r) /
                             std::max(state.vx + param_.tf * state.r, 1e-6)) -
                   delta_r;
  alpha.alpha_rl = std::atan((state.vy - param_.lr * state.r) /
                             std::max(state.vx - param_.tr * state.r, 1e-6));
  alpha.alpha_rr = std::atan((state.vy - param_.lr * state.r) /
                             std::max(state.vx + param_.tr * state.r, 1e-6));

  return alpha;
}

FourWheelSimModel::TireForces FourWheelSimModel::getTireForces(
    const SimState& state, const SlipAngles& alpha, const NormalForces& F_N,
    double D, double B) const {
  TireForces F_tire;

  auto pacejka = [](double a, double B, double C, double D, double E) {
    return D * std::sin(C * std::atan(B * a - E * (B * a - std::atan(B * a))));
  };

  // Lateral forces
  F_tire.Fy_fl = F_N.F_N_fl * param_.Df *
                 pacejka(alpha.alpha_fl, param_.Bf, param_.Cf, 1.0, 0.0);
  F_tire.Fy_fr = F_N.F_N_fr * param_.Df *
                 pacejka(alpha.alpha_fr, param_.Bf, param_.Cf, 1.0, 0.0);
  F_tire.Fy_rl = F_N.F_N_rl * param_.Dr *
                 pacejka(alpha.alpha_rl, param_.Br, param_.Cr, 1.0, 0.0);
  F_tire.Fy_rr = F_N.F_N_rr * param_.Dr *
                 pacejka(alpha.alpha_rr, param_.Br, param_.Cr, 1.0, 0.0);

  // Longitudinal forces
  double F_M = param_.Cm1 * D - param_.Cm2 * D * state.vx;

  // Basic open differential (50/50 split)
  F_tire.Fx_rl = 0.5 * F_M - 0.5 * param_.CBr * B;
  F_tire.Fx_rr = 0.5 * F_M - 0.5 * param_.CBr * B;
  F_tire.Fx_fl = -0.5 * param_.CBf * B;
  F_tire.Fx_fr = -0.5 * param_.CBf * B;

  // Add constant friction
  double F_fric_const = -param_.Cr0;
  F_tire.Fx_rl += 0.25 * F_fric_const;
  F_tire.Fx_rr += 0.25 * F_fric_const;
  F_tire.Fx_fl += 0.25 * F_fric_const;
  F_tire.Fx_fr += 0.25 * F_fric_const;

  return F_tire;
}

SimState FourWheelSimModel::getDerivative(const SimState& s,
                                          const Input& u) const {
  SimState ds;

  AeroForces aero = getAeroForces(s);

  // Approximate accelerations for load transfer
  double F_M_approx = param_.Cm1 * s.D - param_.Cm2 * s.D * s.vx;
  double F_B_approx = -(param_.CBf + param_.CBr) * s.B;
  double ax = (F_M_approx + F_B_approx + aero.F_drag - param_.Cr0) / param_.m;
  double ay = s.vx * s.r;

  NormalForces F_N = getNormalForces(s, aero, ax, ay);
  SlipAngles alpha = getSlipAngles(s, s.delta);
  TireForces F_tire = getTireForces(s, alpha, F_N, s.D, s.B);

  double Fx_f_total = F_tire.Fx_fl + F_tire.Fx_fr;
  double Fx_r_total = F_tire.Fx_rl + F_tire.Fx_rr;
  double Fy_f_total = F_tire.Fy_fl + F_tire.Fy_fr;
  double Fy_r_total = F_tire.Fy_rl + F_tire.Fy_rr;

  // Dynamics
  ds.X = s.vx * std::cos(s.psi) - s.vy * std::sin(s.psi);
  ds.Y = s.vx * std::sin(s.psi) + s.vy * std::cos(s.psi);
  ds.psi = s.r;

  double ax_long = Fx_f_total * std::cos(s.delta) + Fx_r_total + aero.F_drag -
                   Fy_f_total * std::sin(s.delta) + param_.m * s.vy * s.r;

  double ay_lat = Fy_r_total + Fy_f_total * std::cos(s.delta) +
                  Fx_f_total * std::sin(s.delta) - param_.m * s.vx * s.r;

  double domega =
      Fy_f_total * std::cos(s.delta) * param_.lf +
      (F_tire.Fy_fl - F_tire.Fy_fr) * std::sin(s.delta) * param_.tf +
      Fx_f_total * std::sin(s.delta) * param_.lf +
      (F_tire.Fx_fl - F_tire.Fx_fr) * std::cos(s.delta) * param_.tf -
      Fy_r_total * param_.lr + (F_tire.Fx_rl - F_tire.Fx_rr) * param_.tr;

  ds.vx = (1.0 / param_.m) * ax_long;
  ds.vy = (1.0 / param_.m) * ay_lat;
  ds.r = (1.0 / param_.Iz) * domega;

  ds.delta = u.dDelta;
  ds.D = u.dD;
  ds.B = u.dB;

  return ds;
}

State FourWheelSimModel::getMPCCState(double s_approx, double vs) const {
  State x;
  x.X = x_sim_.X;
  x.Y = x_sim_.Y;
  x.phi = x_sim_.psi;
  x.vx = x_sim_.vx;
  x.vy = x_sim_.vy;
  x.r = x_sim_.r;
  x.s = s_approx;
  x.D = x_sim_.D;
  x.B = x_sim_.B;
  x.delta = x_sim_.delta;
  x.vs = vs;
  return x;
}

}  // namespace mpcc
