// Copyright 2019 Alexander Liniger

// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at

//     http://www.apache.org/licenses/LICENSE-2.0

// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#include "Constraints/constraints.h"

#include <iostream>
namespace mpcc {
Constraints::Constraints() {
  std::cout << "default constructor, not everything is initialized properly"
            << std::endl;
}

Constraints::Constraints(double Ts, const PathToJson &path)
    : model_(Ts, path), param_(Param(path.param_path)) {
  tire_con_front_lib_ = (std::make_unique<CppAD::cg::LinuxDynamicLib<double>>(
      path.adcodegen_path + "/cppad_cg_TCF.so"));
  tire_con_rear_lib_ = (std::make_unique<CppAD::cg::LinuxDynamicLib<double>>(
      path.adcodegen_path + "/cppad_cg_TCR.so"));
  tire_con_front_model_ = tire_con_front_lib_->model("TireConFront");
  tire_con_rear_model_ = tire_con_rear_lib_->model("TireConRear");
}

OneDConstraint Constraints::getTrackConstraints(const BoostSplines &track,
                                                const State &x) const {
  // Given arc length s and the track -> compute linearized track constraints
  const double s = x.s;

  // X-Y point of the center line and its derivative
  const Eigen::Vector2d pos_center = track.getPostion(s);
  const Eigen::Vector2d d_center = track.getDerivative(s);

  // Normal vector to the center line (perpendicular to path)
  const Eigen::Vector2d normal = {-d_center(1), d_center(0)};
  const double path_heading = std::atan2(d_center(1), d_center(0));

  // Inner and outer track boundary given left and right width of track
  const double n_left = track.getNLeft(s);    // Positive (usually)
  const double n_right = track.getNRight(s);  // Negative (usually)

  // Alignment of the car with the path
  const double heading_error = std::atan2(std::sin(x.phi - path_heading),
                                          std::cos(x.phi - path_heading));

  const double cos_mu = std::cos(heading_error);
  const double sin_mu = std::sin(heading_error);

  // Projected car size onto the track normal
  const double projected_size =
      param_.car_w * cos_mu + param_.car_l * std::fabs(sin_mu);
  const double d_projected_size =
      -param_.car_w * sin_mu +
      param_.car_l * cos_mu * (heading_error > 0 ? 1.0 : -1.0);

  const Eigen::Vector2d pos_outer =
      pos_center + (n_left - 0.5 * projected_size) * normal;
  const Eigen::Vector2d pos_inner =
      pos_center + (n_right + 0.5 * projected_size) * normal;

  // Compute boundaries for the Taylor expansion
  // dl <= C*x <= du
  const Eigen::Vector2d car_pos = {x.X, x.Y};
  const double current_proj = normal.dot(car_pos);

  const double track_constraint_lower = normal.dot(pos_inner) - current_proj;
  const double track_constraint_upper = normal.dot(pos_outer) - current_proj;

  // Define track Jacobian (direction of the constraint)
  C_i_MPC C_track_constraint = C_i_MPC::Zero();
  C_track_constraint(0, si_index.X) = normal(0);
  C_track_constraint(0, si_index.Y) = normal(1);

  // Add phi dependency to Jacobian based on which bound is closer
  if (std::fabs(track_constraint_lower) < std::fabs(track_constraint_upper)) {
    C_track_constraint(0, si_index.phi) = -0.5 * d_projected_size;
  } else {
    C_track_constraint(0, si_index.phi) = 0.5 * d_projected_size;
  }

  return {C_track_constraint, track_constraint_lower, track_constraint_upper};
}

OneDConstraint Constraints::getTireConstraintRear(const State &x) const {
  // compute linearized slip angle constraints
  // -Inf <= F_comb - F_max <= 0
  const StateVector x_vec = stateToVector(x);
  const std::vector<double> x_std_vec(x_vec.data(),
                                      x_vec.data() + x_vec.size());

  // compute the jacobean of alpha_f
  const C_i_MPC C_tire_constraint_rear =
      Eigen::Map<C_i_MPC>((tire_con_rear_model_->Jacobian(x_std_vec)).data());
  // compute the bounds given the Tylor series expansion
  const double tire_constraint_rear_lower = -INF;
  const double tire_constraint_rear_upper =
      -tire_con_rear_model_->ForwardZero(x_std_vec)[0] -
      C_tire_constraint_rear * x_vec;

  return {C_tire_constraint_rear, tire_constraint_rear_lower,
          tire_constraint_rear_upper};
}

OneDConstraint Constraints::getTireConstraintFront(const State &x) const {
  // compute linearized slip angle constraints
  // -Inf <= F_comb - F_max <= 0
  const StateVector x_vec = stateToVector(x);
  const std::vector<double> x_std_vec(x_vec.data(),
                                      x_vec.data() + x_vec.size());

  // compute the jacobean of alpha_f
  const C_i_MPC C_tire_constraint_front =
      Eigen::Map<C_i_MPC>((tire_con_front_model_->Jacobian(x_std_vec)).data());
  // compute the bounds given the Tylor series expansion
  const double tire_constraint_front_lower = -INF;
  const double tire_constraint_front_upper =
      -tire_con_front_model_->ForwardZero(x_std_vec)[0] -
      C_tire_constraint_front * x_vec;

  return {C_tire_constraint_front, tire_constraint_front_lower,
          tire_constraint_front_upper};
}

ConstrainsMatrix Constraints::getConstraints(const BoostSplines &track,
                                             const State &x,
                                             const Input &u) const {
  // compute all the polytopic state constraints
  // compute the three constraints

  ConstrainsMatrix constrains_matrix;
  const OneDConstraint track_constraints = getTrackConstraints(track, x);
  const OneDConstraint tire_constraints_rear = getTireConstraintRear(x);
  const OneDConstraint tire_constraints_front = getTireConstraintFront(x);

  C_MPC C_constrains_matrix;
  d_MPC dl_constrains_matrix;
  d_MPC du_constrains_matrix;

  C_constrains_matrix.row(si_index.con_track) = track_constraints.C_i;
  dl_constrains_matrix(si_index.con_track) = track_constraints.dl_i;
  du_constrains_matrix(si_index.con_track) = track_constraints.du_i;

  C_constrains_matrix.row(si_index.con_tire_r) = tire_constraints_rear.C_i;
  dl_constrains_matrix(si_index.con_tire_r) = tire_constraints_rear.dl_i;
  du_constrains_matrix(si_index.con_tire_r) = tire_constraints_rear.du_i;

  C_constrains_matrix.row(si_index.con_tire_f) = tire_constraints_front.C_i;
  dl_constrains_matrix(si_index.con_tire_f) = tire_constraints_front.dl_i;
  du_constrains_matrix(si_index.con_tire_f) = tire_constraints_front.du_i;

  return {C_constrains_matrix, D_MPC::Zero(), dl_constrains_matrix,
          du_constrains_matrix};
}
}  // namespace mpcc
