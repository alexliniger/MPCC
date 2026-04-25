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

#ifndef MPCC_CONSTRAINTS_H
#define MPCC_CONSTRAINTS_H

#include "config.h"
#include "Spline/arc_length_spline.h"
#include "Model/model.h"
namespace mpcc{
struct ConstrainsMatrix {
    // dl <= C xk + D uk <= du
    C_MPC C;    //polytopic state constraints
    D_MPC D;    //polytopic input constraints
    d_MPC dl;   //lower bounds
    d_MPC du;   //upper bounds
};

struct OneDConstraint {
    const C_i_MPC C_i;
    const double dl_i;
    const double du_i;
};

class Constraints {
public:
    ConstrainsMatrix getConstraints(const ArcLengthSpline &track,const State &x,const Input &u) const;

    Constraints();
    Constraints(double Ts,const PathToJson &path);
private:
    OneDConstraint getTrackConstraints(const ArcLengthSpline &track,const State &x) const;

    OneDConstraint getTireConstraintRear(const State &x) const;
    C_i_MPC getTireConstraintRearJac(const State &x) const;

    OneDConstraint getAlphaConstraintFront(const State &x) const;
    C_i_MPC getAlphaConstraintFrontJac(const State &x) const;

    Model model_;
    Param param_;
};
}

#endif //MPCC_CONSTRAINTS_H
