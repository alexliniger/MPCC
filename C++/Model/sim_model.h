// Copyright 2019 Alexander Liniger
// Licensed under the Apache License, Version 2.0

#ifndef MPCC_SIM_MODEL_H
#define MPCC_SIM_MODEL_H

#include <string>
#include <vector>

#include "../Config/params.h"
#include "../Config/types.h"

namespace mpcc {

// Interface for stateful simulation models
class ISimModel {
 public:
  virtual ~ISimModel() = default;

  // Initialize the internal complex state from a simple MPC state
  virtual void initialize(const State& x) = 0;

  // Step the internal complex state forward
  virtual void step(const Input& input, double dt) = 0;

  // Extract the rigid body position for track projection
  virtual double getX() const = 0;
  virtual double getY() const = 0;

  // Extract the simple MPC state from the current internal complex state
  virtual State getMPCCState(double s_approx, double vs) const = 0;
};

// -----------------------------------------------------------------------------
// Bicycle Model Implementation
// -----------------------------------------------------------------------------

struct SimState {
  double X, Y, psi;
  double vx, vy, r;
  // Control states integrated from rates
  double delta;
  double D;
  double B;
};

class BicycleSimModel : public ISimModel {
 public:
  BicycleSimModel(double Ts, const std::string& param_path);

  void initialize(const State& x) override;
  void step(const Input& input, double dt) override;

  double getX() const override {
    return x_sim_.X;
  }
  double getY() const override {
    return x_sim_.Y;
  }

  State getMPCCState(double s_approx, double vs) const override;

 private:
  struct NormalForces {
    double F_N_front;
    double F_N_rear;
  };

  struct AeroForces {
    double F_downforce;
    double F_drag;
  };

  struct SlipAngles {
    double alpha_f;
    double alpha_r;
  };

  struct TireForces {
    double Fy_f;
    double Fy_r;
    double Fx_f;
    double Fx_r;
  };

  Param param_;
  double Ts_;
  SimState x_sim_;  // Internal state

  AeroForces getAeroForces(const SimState& state) const;
  NormalForces getNormalForces(const SimState& state,
                               const AeroForces& aero) const;
  SlipAngles getSlipAngles(const SimState& state) const;
  TireForces getTireForces(const SimState& state, const SlipAngles& alpha,
                           const NormalForces& F_N) const;

  SimState getDerivative(const SimState& s, const Input& u) const;
};

}  // namespace mpcc

#endif  // MPCC_SIM_MODEL_H
