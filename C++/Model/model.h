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

#ifndef MPCC_MODEL_H
#define MPCC_MODEL_H

#include "config.h"
#include "types.h"
#include "Params/params.h"
#include <cppad/cg.hpp>

namespace mpcc{
//Return
struct LinModelMatrix {
    A_MPC A;
    B_MPC B;
    g_MPC g;
};

struct TireForces {
    const double F_y;
    const double F_x;
};

struct NormalForces {
    const double F_N_front;
    const double F_N_rear;
};

struct TireForcesDerivatives{
    const double dF_y_vx;
    const double dF_y_vy;
    const double dF_y_r;
    const double dF_y_D;
    const double dF_y_delta;

    const double dF_x_vx;
    const double dF_x_vy;
    const double dF_x_r;
    const double dF_x_D;
    const double dF_x_delta;
};

struct FrictionForceDerivatives {
    const double dF_f_vx;
    const double dF_f_vy;
    const double dF_f_r;
    const double dF_f_D;
    const double dF_f_delta;
};

class Model {
public:
    StateVector getF(const State &x,const Input &u) const;

    LinModelMatrix getLinModel(const State &x, const Input &u,const State &x_next) const;

    Model();
    Model(double Ts,const PathToJson &path);
private:
    LinModelMatrix discretizeModel(const State &x, const Input &u,const State &x_next) const;

    std::unique_ptr<CppAD::cg::LinuxDynamicLib<double>> RK4_lib_;
    std::unique_ptr<CppAD::cg::GenericModel<double>> RK4_model_;

    std::unique_ptr<CppAD::cg::LinuxDynamicLib<double>> f_dyn_lib_;
    std::unique_ptr<CppAD::cg::GenericModel<double>> f_dyn_model_;

    Param param_;
    const double Ts_;
};
}
#endif //MPCC_MODEL_H
