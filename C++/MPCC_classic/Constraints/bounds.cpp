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

#include "bounds.h"
namespace mpcc{
Bounds::Bounds()
{
    std::cout << "default constructor, not everything is initialized properly" << std::endl;
}

Bounds::Bounds(BoundsParam bounds_param) 
{
    l_bounds_x_(0) = bounds_param.lower_state_bounds.X_l;
    l_bounds_x_(1) = bounds_param.lower_state_bounds.Y_l;
    l_bounds_x_(2) = bounds_param.lower_state_bounds.phi_l;
    l_bounds_x_(3) = bounds_param.lower_state_bounds.vx_l;
    l_bounds_x_(4) = bounds_param.lower_state_bounds.vy_l;
    l_bounds_x_(5) = bounds_param.lower_state_bounds.r_l;
    l_bounds_x_(6) = bounds_param.lower_state_bounds.s_l;
    l_bounds_x_(7) = bounds_param.lower_state_bounds.D_l;
    l_bounds_x_(8) = bounds_param.lower_state_bounds.delta_l;
    l_bounds_x_(9) = bounds_param.lower_state_bounds.vs_l;

    u_bounds_x_(0) = bounds_param.upper_state_bounds.X_u;
    u_bounds_x_(1) = bounds_param.upper_state_bounds.Y_u;
    u_bounds_x_(2) = bounds_param.upper_state_bounds.phi_u;
    u_bounds_x_(3) = bounds_param.upper_state_bounds.vx_u;
    u_bounds_x_(4) = bounds_param.upper_state_bounds.vy_u;
    u_bounds_x_(5) = bounds_param.upper_state_bounds.r_u;
    u_bounds_x_(6) = bounds_param.upper_state_bounds.s_u;
    u_bounds_x_(7) = bounds_param.upper_state_bounds.D_u;
    u_bounds_x_(8) = bounds_param.upper_state_bounds.delta_u;
    u_bounds_x_(9) = bounds_param.upper_state_bounds.vs_u;

    l_bounds_u_(0) = bounds_param.lower_input_bounds.dD_l;
    l_bounds_u_(1) = bounds_param.lower_input_bounds.dDelta_l;
    l_bounds_u_(2) = bounds_param.lower_input_bounds.dVs_l;

    u_bounds_u_(0) = bounds_param.upper_input_bounds.dD_u;
    u_bounds_u_(1) = bounds_param.upper_input_bounds.dDelta_u;
    u_bounds_u_(2) = bounds_param.upper_input_bounds.dVs_u;

    l_bounds_s_ = Bounds_s::Zero();
    u_bounds_s_ = Bounds_s::Zero();

    std::cout << "bounds initialized" << std::endl;
}

Bounds_x Bounds::getBoundsLX() const
{
    return  l_bounds_x_;
}

Bounds_x Bounds::getBoundsUX() const
{
    return  u_bounds_x_;
}

Bounds_u Bounds::getBoundsLU() const
{
    return  l_bounds_u_;
}

Bounds_u Bounds::getBoundsUU() const
{
    return  u_bounds_u_;
}

Bounds_s Bounds::getBoundsLS() const
{
    return  l_bounds_s_;
}

Bounds_s Bounds::getBoundsUS() const{
    return  u_bounds_s_;
}
}