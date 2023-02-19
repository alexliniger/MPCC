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

#ifndef MPCC_CONFIG_H
#define MPCC_CONFIG_H

#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>

namespace mpcc{

// #define MAX(a,b) (a < b) ? b : a

#define NX 11
#define NU 4

#define NB 15 //max number of bounds
#define NPC 3 //number of polytopic constraints
#define NS 3

static constexpr int N = 100;
static constexpr double INF = 1E5;
static constexpr int N_SPLINE = 5000;


struct StateInputIndex{
    int X = 0;
    int Y = 1;
    int phi = 2;
    int vx = 3;
    int vy = 4;
    int r = 5;
    int s = 6;
    int D = 7;
    int B = 8;
    int delta = 9;
    int vs = 10;

    int dD = 0;
    int dB = 1;
    int dDelta = 2;
    int dVs = 3;

    int con_track = 0;
    int con_tire_f = 1;
    int con_tire_r = 2;
};

static const StateInputIndex si_index;

}
#endif //MPCC_CONFIG_H
