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

}  // namespace mpcc
