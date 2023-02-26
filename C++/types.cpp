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

#include "types.h"
namespace mpcc{

StateVector stateToVector(const State &x)
{
    StateVector xk;
    xk(si_index.X) = x.X;
    xk(si_index.Y) = x.Y;
    xk(si_index.phi) = x.phi;
    xk(si_index.vx) = x.vx;
    xk(si_index.vy) = x.vy;
    xk(si_index.r) = x.r;
    xk(si_index.s) = x.s;
    xk(si_index.D) = x.D;
    xk(si_index.B) = x.B;
    xk(si_index.delta) = x.delta;
    xk(si_index.vs) = x.vs;
    return xk;
}

InputVector inputToVector(const Input &u)
{
    InputVector uk = {u.dD,u.dB,u.dDelta,u.dVs};
    return uk;
}

State vectorToState(const StateVector &xk)
{
    State x;
    x.X     = xk(si_index.X);
    x.Y     = xk(si_index.Y);
    x.phi   = xk(si_index.phi);
    x.vx    = xk(si_index.vx);
    x.vy    = xk(si_index.vy);
    x.r     = xk(si_index.r);
    x.s     = xk(si_index.s);
    x.D     = xk(si_index.D);
    x.B     = xk(si_index.B);
    x.delta = xk(si_index.delta);
    x.vs    = xk(si_index.vs);

    return x;
}

Input vectorToInput(const InputVector &uk)
{
    Input u;
    u.dD     = uk(si_index.dD);
    u.dB     = uk(si_index.dB);
    u.dDelta = uk(si_index.dDelta);
    u.dVs    = uk(si_index.dVs);

    return u;
}

State arrayToState(double *xk)
{
    State x;
    x.X     = xk[si_index.X];
    x.Y     = xk[si_index.Y];
    x.phi   = xk[si_index.phi];
    x.vx    = xk[si_index.vx];
    x.vy    = xk[si_index.vy];
    x.r     = xk[si_index.r];
    x.s     = xk[si_index.s];
    x.D     = xk[si_index.D];
    x.B     = xk[si_index.B];
    x.delta = xk[si_index.delta];
    x.vs    = xk[si_index.vs];

    return x;
}

Input arrayToInput(double *uk)
{
    Input u;
    u.dD     = uk[si_index.dD];
    u.dB     = uk[si_index.dB];
    u.dDelta = uk[si_index.dDelta];
    u.dVs    = uk[si_index.dVs];

    return u;
}

std::vector<double> stateInputToVector(const State x,const Input u)
{

    StateVector xv = stateToVector(x);
    InputVector uv = inputToVector(u);

    Eigen::Vector<double,NX+NU> z;
    z << xv,uv;
    std::vector<double> zv(z.data(), z.data() + z.size());

    return zv;
}

}