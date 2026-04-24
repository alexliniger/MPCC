//
// Created by alexliniger on 24.07.20.
//

#include "ADCodeGen/ad_dynamics.h"

#include <iostream>
#include <memory>
#include <string>
#include <vector>
namespace mpcc {

using CppAD::ADFun;
using CppAD::Independent;
using CppAD::cg::DynamicModelLibraryProcessor;
using CppAD::cg::GccCompiler;
using CppAD::cg::ModelCSourceGen;
using CppAD::cg::ModelLibraryCSourceGen;
using CppAD::cg::SaveFilesModelLibraryProcessor;

template <typename T>
using DynamicLib = CppAD::cg::DynamicLib<T>;

ADDynamics::ADDynamics() : Ts_(1.0) {
  std::cout << "default constructor, not everything is initialized properly"
            << std::endl;
}
ADDynamics::ADDynamics(double Ts, const std::string& path,
                       const std::string& outputPath)
    : Ts_(Ts), param_(path), outputPath_(outputPath) {
  if (outputPath_.back() != '/' && outputPath_.back() != '\\') {
    outputPath_ += "/";
  }
}
void ADDynamics::genLibraryIntegrator(IntegratorType type, int n_steps) {
  // independent variable vector
  std::vector<ADCG> z(NX + NU);
  for (int i = 0; i < NX + NU; i++) z[i] = 0.0;
  z[3] = 10.;
  Independent(z);

  std::vector<ADCG> x(NX);
  std::vector<ADCG> u(NU);
  for (int i = 0; i < NX; i++) x[i] = z[i];
  for (int i = 0; i < NU; i++) u[i] = z[i + NX];

  // Integrate
  std::vector<ADCG> x_plus(NX);
  x_plus = Integrate(x, u, type, n_steps);
  ADFun<CGD> fun(z, x_plus);

  std::string name = (type == IntegratorType::RK4) ? "RK4" : "Euler";
  std::string lib_name = "cppad_cg_" + name;

  // generates source code
  ModelCSourceGen<double> cgen(fun, name);
  cgen.setCreateJacobian(true);
  ModelLibraryCSourceGen<double> libcgen(cgen);

  // compile source code
  DynamicModelLibraryProcessor<double> p(libcgen, outputPath_ + lib_name);

  GccCompiler<double> compiler;
  std::unique_ptr<DynamicLib<double>> dynamicLib =
      p.createDynamicLibrary(compiler);

  // save to files (not really required)
  SaveFilesModelLibraryProcessor<double> p2(libcgen);
  p2.saveSourcesTo(outputPath_ + "sources");
}

void ADDynamics::genLibraryGetF() {
  // independent variable vector
  std::vector<ADCG> z(NX + NU);
  for (int i = 0; i < NX + NU; i++) z[i] = 0.0;
  z[3] = 10.;
  Independent(z);
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // RK4
  // dependent variable vector
  std::vector<ADCG> x_dot(NX);
  // the model equation
  x_dot = f_dyn(z);
  ADFun<CGD> fun(z, x_dot);

  // generates source code
  ModelCSourceGen<double> cgen(fun, "f_dyn");
  cgen.setCreateJacobian(true);
  ModelLibraryCSourceGen<double> libcgen(cgen);

  // compile source code
  DynamicModelLibraryProcessor<double> p(libcgen,
                                         outputPath_ + "cppad_cg_f_dyn");

  GccCompiler<double> compiler;
  std::unique_ptr<DynamicLib<double>> dynamicLibRK4 =
      p.createDynamicLibrary(compiler);

  // save to files (not really required)
  SaveFilesModelLibraryProcessor<double> p2(libcgen);
  p2.saveSourcesTo(outputPath_ + "sources");
}

void ADDynamics::genLibraryTireFront() {
  // Front Combined Force
  std::vector<ADCG> x(NX);
  for (int i = 0; i < NX; i++) x[i] = 0.0;
  x[3] = 10.;
  Independent(x);
  // dependent variable vector
  std::vector<ADCG> F_comb_f(1);
  // the model equation
  F_comb_f = TireConFront(x);
  ADFun<CGD> fun(x, F_comb_f);

  // generates source code
  ModelCSourceGen<double> cgen(fun, "TireConFront");
  cgen.setCreateJacobian(true);
  ModelLibraryCSourceGen<double> libcgen(cgen);

  // compile source code
  DynamicModelLibraryProcessor<double> p(libcgen, outputPath_ + "cppad_cg_TCF");

  GccCompiler<double> compiler;
  std::unique_ptr<DynamicLib<double>> dynamicLibTireConFront =
      p.createDynamicLibrary(compiler);

  // save to files (not really required)
  SaveFilesModelLibraryProcessor<double> p2(libcgen);
  p2.saveSourcesTo(outputPath_ + "sources");
}

void ADDynamics::genLibraryTireRear() {
  // Front Combined Force
  std::vector<ADCG> x(NX);
  for (int i = 0; i < NX; i++) x[i] = 0.0;
  x[3] = 10.;
  Independent(x);
  // dependent variable vector
  std::vector<ADCG> F_comb_r(1);
  // the model equation
  F_comb_r = TireConRear(x);
  ADFun<CGD> fun(x, F_comb_r);

  // generates source code
  ModelCSourceGen<double> cgen(fun, "TireConRear");
  cgen.setCreateJacobian(true);
  ModelLibraryCSourceGen<double> libcgen(cgen);

  // compile source code
  DynamicModelLibraryProcessor<double> p(libcgen, outputPath_ + "cppad_cg_TCR");

  GccCompiler<double> compiler;
  std::unique_ptr<DynamicLib<double>> dynamicLibTireConRear =
      p.createDynamicLibrary(compiler);

  // save to files (not really required)
  SaveFilesModelLibraryProcessor<double> p2(libcgen);
  p2.saveSourcesTo(outputPath_ + "sources");
}

std::vector<ADCG> ADDynamics::scalerMult(std::vector<ADCG> x, double a) {
  std::vector<ADCG> y(NX);
  for (int i = 0; i < NX; i++) y[i] = a * x[i];
  return y;
}

std::vector<ADCG> ADDynamics::vectorAdd(std::vector<ADCG> x1,
                                        std::vector<ADCG> x2) {
  std::vector<ADCG> y(NX);
  for (int i = 0; i < NX; i++) y[i] = x1[i] + x2[i];
  return y;
}

ADCG ADDynamics::getState(std::vector<ADCG> x, int index) {
  return x[index];
}

ADCG ADDynamics::getInput(std::vector<ADCG> u, int index) {
  return u[index];
}

ADCG ADDynamics::getSlipAngleFront(std::vector<ADCG> x) {
  ADCG vx = getState(x, si_index.vx);
  // compute slip angels given current state
  return atan((getState(x, si_index.vy) + getState(x, si_index.r) * param_.lf) /
              vx) -
         getState(x, si_index.delta);
}

ADCG ADDynamics::getSlipAngleRear(std::vector<ADCG> x) {
  // compute slip angels given current state
  ADCG vx = getState(x, si_index.vx);
  return atan((getState(x, si_index.vy) - getState(x, si_index.r) * param_.lr) /
              vx);
}

TireForces ADDynamics::getForceFront(std::vector<ADCG> x) {
  ADCG alpha_f = getSlipAngleFront(x);
  NormalForces f_normal = getForceNormalDyn(x);

  ADCG F_y = f_normal.F_N_front * param_.Df *
             sin(param_.Cf * atan(param_.Bf * alpha_f));
  ADCG F_x = -param_.CBf * getState(x, si_index.B) - param_.Cr0 * 0.5;

  return {F_x, F_y};
}

TireForces ADDynamics::getForceRear(std::vector<ADCG> x) {
  ADCG alpha_r = getSlipAngleRear(x);
  NormalForces f_normal = getForceNormalDyn(x);

  ADCG F_y = f_normal.F_N_rear * param_.Dr *
             sin(param_.Cr * atan(param_.Br * alpha_r));
  ADCG F_x = param_.Cm1 * getState(x, si_index.D) -
             param_.Cm2 * getState(x, si_index.D) * getState(x, si_index.vx) -
             param_.CBr * getState(x, si_index.B) - param_.Cr0 * 0.5;

  return {F_x, F_y};
}

ADCG ADDynamics::getForceFriction(std::vector<ADCG> x) {
  return -0.5 * param_.rho * param_.S * param_.Cr2 * getState(x, si_index.vx) *
         getState(x, si_index.vx);
}

NormalForces ADDynamics::getForceNormalStatic(void) {
  ADCG F_N_front =
      (ADCG)(param_.lr / (param_.lf + param_.lr) * param_.m * param_.g);
  ADCG F_N_rear =
      (ADCG)(param_.lf / (param_.lf + param_.lr) * param_.m * param_.g);
  return {F_N_front, F_N_rear};
}

NormalForces ADDynamics::getForceNormalDyn(std::vector<ADCG> x) {
  // including aero
  ADCG vx = getState(x, si_index.vx);

  ADCG F_N_front = param_.aero_split_front * 0.5 * param_.rho * param_.S *
                       param_.Cl * vx * vx +
                   param_.lr / (param_.lf + param_.lr) * param_.m * param_.g;
  ADCG F_N_rear = (1.0 - param_.aero_split_front) * 0.5 * param_.rho *
                      param_.S * param_.Cl * vx * vx +
                  param_.lf / (param_.lf + param_.lr) * param_.m * param_.g;
  return {F_N_front, F_N_rear};
}

std::vector<ADCG> ADDynamics::dx(std::vector<ADCG> x, std::vector<ADCG> u) {
  ADCG phi = getState(x, si_index.phi);
  ADCG vx = getState(x, si_index.vx);
  ADCG vy = getState(x, si_index.vy);
  ADCG r = getState(x, si_index.r);
  ADCG D = getState(x, si_index.D);
  ADCG delta = getState(x, si_index.delta);
  ADCG vs = getState(x, si_index.vs);

  ADCG dD = getInput(u, si_index.dD);
  ADCG dB = getInput(u, si_index.dB);
  ADCG dDelta = getInput(u, si_index.dDelta);
  ADCG dVs = getInput(u, si_index.dVs);

  TireForces tire_forces_front = getForceFront(x);
  TireForces tire_forces_rear = getForceRear(x);
  ADCG friction_force = getForceFriction(x);

  std::vector<ADCG> dx(NX);
  dx[0] = vx * cos(phi) - vy * sin(phi);
  dx[1] = vy * cos(phi) + vx * sin(phi);
  dx[2] = r;
  dx[3] = 1.0 / param_.m *
          (tire_forces_rear.F_x + friction_force -
           tire_forces_front.F_y * sin(delta) +
           tire_forces_front.F_x * cos(delta) + param_.m * vy * r);
  dx[4] = 1.0 / param_.m *
          (tire_forces_rear.F_y + tire_forces_front.F_x * sin(delta) +
           tire_forces_front.F_y * cos(delta) - param_.m * vx * r);
  dx[5] = 1.0 / param_.Iz *
          (-tire_forces_rear.F_y * param_.lr +
           (tire_forces_front.F_x * sin(delta) +
            tire_forces_front.F_y * cos(delta)) *
               param_.lf);
  dx[6] = vs;
  dx[7] = dD;
  dx[8] = dB;
  dx[9] = dDelta;
  dx[10] = dVs;

  return dx;
}

std::vector<ADCG> ADDynamics::RK4(std::vector<ADCG> x, std::vector<ADCG> u,
                                  double dt) {
  std::vector<ADCG> k1 = dx(x, u);
  std::vector<ADCG> k2 = dx(vectorAdd(x, scalerMult(k1, dt * 0.5)), u);
  std::vector<ADCG> k3 = dx(vectorAdd(x, scalerMult(k2, dt * 0.5)), u);
  std::vector<ADCG> k4 = dx(vectorAdd(x, scalerMult(k3, dt)), u);

  std::vector<ADCG> y(NX);
  for (int i = 0; i < NX; i++)
    y[i] = x[i] + dt * (k1[i] / 6.0 + k2[i] / 3.0 + k3[i] / 3.0 + k4[i] / 6.0);

  return y;
}

std::vector<ADCG> ADDynamics::ForwardEuler(std::vector<ADCG> x,
                                           std::vector<ADCG> u, double dt) {
  std::vector<ADCG> k1 = dx(x, u);
  return vectorAdd(x, scalerMult(k1, dt));
}

std::vector<ADCG> ADDynamics::Integrate(std::vector<ADCG> x,
                                        std::vector<ADCG> u,
                                        IntegratorType type, int n_steps) {
  std::vector<ADCG> x_next = x;
  double dt = Ts_ / static_cast<double>(n_steps);
  for (int i = 0; i < n_steps; i++) {
    if (type == IntegratorType::RK4) {
      x_next = RK4(x_next, u, dt);
    } else {
      x_next = ForwardEuler(x_next, u, dt);
    }
  }
  return x_next;
}
std::vector<ADCG> ADDynamics::TireConFront(std::vector<ADCG> x) {
  TireForces tire_forces_front = getForceFront(x);
  NormalForces f_normal = getForceNormalStatic();
  NormalForces f_normal_dyn = getForceNormalDyn(x);

  ADCG comb_force = ((param_.e_long * param_.e_long * tire_forces_front.F_x *
                      tire_forces_front.F_x) +
                     (tire_forces_front.F_y * tire_forces_front.F_y) -
                     (param_.e_eps * param_.e_eps * f_normal_dyn.F_N_front *
                      f_normal_dyn.F_N_front * param_.Df * param_.Df)) /
                    (f_normal.F_N_front * f_normal.F_N_front);

  return {comb_force};
}

std::vector<ADCG> ADDynamics::TireConRear(std::vector<ADCG> x) {
  TireForces tire_forces_rear = getForceRear(x);
  NormalForces f_normal = getForceNormalStatic();
  NormalForces f_normal_dyn = getForceNormalDyn(x);

  ADCG comb_force = ((param_.e_long * param_.e_long * tire_forces_rear.F_x *
                      tire_forces_rear.F_x) +
                     (tire_forces_rear.F_y * tire_forces_rear.F_y) -
                     (param_.e_eps * param_.e_eps * f_normal_dyn.F_N_rear *
                      f_normal_dyn.F_N_rear * param_.Dr * param_.Dr)) /
                    (f_normal.F_N_rear * f_normal.F_N_rear);

  return {comb_force};
}
std::vector<ADCG> ADDynamics::f_dyn(std::vector<ADCG> z) {
  std::vector<ADCG> state(NX);
  std::vector<ADCG> input(NU);
  for (int i = 0; i < NX; i++) {
    state[i] = z[i];
  }
  for (int i = 0; i < NU; i++) {
    input[i] = z[i + NX];
  }

  return dx(state, input);
}

}  // namespace mpcc
