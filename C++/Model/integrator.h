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

#ifndef MPCC_INTEGRATOR_H
#define MPCC_INTEGRATOR_H

#include "config.h"
#include "model.h"
#include "types.h"

namespace mpcc{
class Integrator {
public:
    State RK4(const State &x, const Input &u,double ts) const;
    State EF(const State &x, const Input &u,double ts) const;
    State simTimeStep(const State &x, const Input &u,double ts) const;

    Integrator();
    Integrator(double Ts, const PathToJson &path);

private:
    const double fine_time_step_ = 0.001;

    Model model_;
};
}
#endif //MPCC_INTEGRATOR_H
