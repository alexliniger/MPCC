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

#ifndef MPCC_TYPES_H
#define MPCC_TYPES_H

#include "config.h"
namespace mpcc{
struct State{
    double X;
    double Y;
    double phi;
    double vx;
    double vy;
    double r;
    double s;
    double D;
    double delta;
    double vs;

    void setZero()
    {
        X = 0.0;
        Y = 0.0;
        phi = 0.0;
        vx = 0.0;
        vy = 0.0;
        r = 0.0;
        s = 0.0;
        D = 0.0;
        delta = 0.0;
        vs = 0.0;
    }

    void unwrap(double track_length)
    {
        if (phi > M_PI)
            phi -= 2.0 * M_PI;
        if (phi < -M_PI)
            phi += 2.0 * M_PI;

        if (s > track_length)
            s -= track_length;
        if (s < 0)
            s += track_length;
    }

    void vxNonZero(double vx_zero)
    {
        if(vx < vx_zero){
            vx = vx_zero;
            vy = 0.0;
            r = 0.0;
            delta = 0.0;
        }
    }
};

struct Input{
    double dD;
    double dDelta;
    double dVs;

    void setZero()
    {
        dD = 0.0;
        dDelta = 0.0;
        dVs = 0.0;
    }
};

struct PathToJson{
    const std::string param_path;
    const std::string cost_path;
    const std::string bounds_path;
    const std::string track_path;
    const std::string normalization_path;
};

typedef Eigen::Matrix<double,NX,1> StateVector;
typedef Eigen::Matrix<double,NU,1> InputVector;

typedef Eigen::Matrix<double,NX,NX> A_MPC;
typedef Eigen::Matrix<double,NX,NU> B_MPC;
typedef Eigen::Matrix<double,NX,1> g_MPC;

typedef Eigen::Matrix<double,NX,NX> Q_MPC;
typedef Eigen::Matrix<double,NU,NU> R_MPC;
typedef Eigen::Matrix<double,NX,NU> S_MPC;

typedef Eigen::Matrix<double,NX,1> q_MPC;
typedef Eigen::Matrix<double,NU,1> r_MPC;

typedef Eigen::Matrix<double,NPC,NX> C_MPC;
typedef Eigen::Matrix<double,1,NX> C_i_MPC;
typedef Eigen::Matrix<double,NPC,NU> D_MPC;
typedef Eigen::Matrix<double,NPC,1> d_MPC;

typedef Eigen::Matrix<double,NS,NS> Z_MPC;
typedef Eigen::Matrix<double,NS,1> z_MPC;

typedef Eigen::Matrix<double,NX,NX> TX_MPC;
typedef Eigen::Matrix<double,NU,NU> TU_MPC;
typedef Eigen::Matrix<double,NS,NS> TS_MPC;

typedef Eigen::Matrix<double,NX,1> Bounds_x;
typedef Eigen::Matrix<double,NU,1> Bounds_u;
typedef Eigen::Matrix<double,NS,1> Bounds_s;

StateVector stateToVector(const State &x);
InputVector inputToVector(const Input &u);

State vectorToState(const StateVector &xk);
Input vectorToInput(const InputVector &uk);

State arrayToState(double *xk);
Input arrayToInput(double *uk);
}
#endif //MPCC_TYPES_H
