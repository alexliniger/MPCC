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

#include "integrator.h"
namespace mpcc{
Integrator::Integrator()
{
    std::cout << "default constructor, not everything is initialized properly" << std::endl;
}

Integrator::Integrator(double Ts,const PathToJson &path)
:model_(Ts,path)
{
}

State Integrator::RK4(const State &x, const Input &u,const double ts) const
{
    // 4th order Runge Kutta (RK4) implementation
    // 4 evaluation points of continuous dynamics
    const StateVector x_vec = stateToVector(x);
    const InputVector u_vec = inputToVector(u);
    // evaluating the 4 points
    const StateVector k1 = model_.getF(vectorToState(x_vec),u);
    const StateVector k2 = model_.getF(vectorToState(x_vec+ts/2.*k1),u);
    const StateVector k3 = model_.getF(vectorToState(x_vec+ts/2.*k2),u);
    const StateVector k4 = model_.getF(vectorToState(x_vec+ts*k3),u);
    // combining to give output
    const StateVector x_next = x_vec + ts*(k1/6.+k2/3.+k3/3.+k4/6.);
    return vectorToState(x_next);
}

State Integrator::EF(const State &x, const Input &u,const double ts) const
{
    // Euler Forward integration
    const StateVector x_vec = stateToVector(x);
    const StateVector f = model_.getF(x,u); // evaluation continuous dynmaics
    // compute next time step
    const StateVector x_next = x_vec + ts*f;
    return vectorToState(x_next);
}

State Integrator::simTimeStep(const State &x, const Input &u,const double ts) const
{
    // integrate time step
    State x_next = x;
    const int integration_steps = (int)(ts/fine_time_step_);
    if(ts/fine_time_step_ != integration_steps)
    {
        std::cout << "Warning" << std::endl;
    }
    for(int i = 0;i<integration_steps;i++)
        x_next = RK4(x_next,u,fine_time_step_);

    return x_next;
}
}