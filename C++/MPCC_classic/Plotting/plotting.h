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

#ifndef MPCC_PLOTTING_H
#define MPCC_PLOTTING_H

#include "config.h"
#include "types.h"
#include "Params/params.h"
#include "Params/track.h"
#include <matplotlibcpp.h>
#include <vector>
#include <MPC/mpc.h>
#include <Model/model.h>

namespace plt = matplotlibcpp;

namespace mpcc {
class Plotting {
public:
    void plotRun(const std::list<MPCReturn> &log, const TrackPos &track_xy) const;
    void plotSim(const std::list<MPCReturn> &log, const TrackPos &track_xy) const;

    Plotting(double Ts, PathToJson path);

private:
    void plotBox(const State &x0) const;

    Model model_;
    Param param_;
};
}

#endif //MPCC_PLOTTING_H
