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

class Model {
public:
    double getSlipAngleFront(const State &x) const;
    double getSlipAngleRear(const State &x) const;

    TireForces getForceFront(const State &x) const;
    TireForces getForceRear(const State &x) const;

    TireForcesDerivatives getForceFrontDerivatives(const State &x) const;
    TireForcesDerivatives getForceRearDerivatives(const State &x) const;

    StateVector getF(const State &x,const Input &u) const;

    LinModelMatrix getLinModel(const State &x, const Input &u) const;
private:
    LinModelMatrix getModelJacobian(const State &x, const Input &u) const;
    LinModelMatrix discretizeModel(const LinModelMatrix &lin_model_c) const;
};
}
#endif //MPCC_MODEL_H
