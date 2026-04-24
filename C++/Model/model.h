// Copyright 2019 Alexander Liniger
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////

#ifndef MPCC_MODEL_H
#define MPCC_MODEL_H

#include <memory>

#include <cppad/cg.hpp>

#include "Config/config.h"
#include "Config/params.h"
#include "Config/types.h"

namespace mpcc {
// Return
struct LinModelMatrix {
  A_MPC A;
  B_MPC B;
  g_MPC g;
};

class Model {
 public:
  StateVector getF(const State &x, const Input &u) const;

  LinModelMatrix getLinModel(const State &x, const Input &u,
                             const State &x_next) const;

  Model();
  Model(double Ts, const PathToJson &path);

 private:
  LinModelMatrix discretizeModel(const State &x, const Input &u,
                                 const State &x_next) const;

  std::unique_ptr<CppAD::cg::LinuxDynamicLib<double>> integrator_lib_;
  std::unique_ptr<CppAD::cg::GenericModel<double>> integrator_model_;

  std::unique_ptr<CppAD::cg::LinuxDynamicLib<double>> f_dyn_lib_;
  std::unique_ptr<CppAD::cg::GenericModel<double>> f_dyn_model_;

  Param param_;
  const double Ts_;
};
}  // namespace mpcc
#endif  // MPCC_MODEL_H
