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
    double getSlipAngleFront(const State &x) const;
    double getSlipAngleRear(const State &x) const;

    TireForces getForceFront(const State &x) const;
    TireForces getForceRear(const State &x) const;
    double getForceFriction(const State &x) const;
    NormalForces getForceNormal(const State &x) const;

    TireForcesDerivatives getForceFrontDerivatives(const State &x) const;
    TireForcesDerivatives getForceRearDerivatives(const State &x) const;
    FrictionForceDerivatives getForceFrictionDerivatives(const State &x) const;

    StateVector getF(const State &x,const Input &u) const;

    LinModelMatrix getLinModel(const State &x, const Input &u) const;

    Model();
    Model(double Ts,const PathToJson &path);
private:
    LinModelMatrix getModelJacobian(const State &x, const Input &u) const;
    LinModelMatrix discretizeModel(const LinModelMatrix &lin_model_c) const;

    Param param_;
    const double Ts_;
};
}
#endif //MPCC_MODEL_H
