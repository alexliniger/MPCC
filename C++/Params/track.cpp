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

#include "track.h"
namespace mpcc{
TrackPos Track::getTrack()
{
    /////////////////////////////////////////////////////
    // Loading Model and Constraint Parameters //////////
    /////////////////////////////////////////////////////
    std::ifstream iTrack("Params/track.json");
    json jsonTrack;
    iTrack >> jsonTrack;
    // Model Parameters
    std::vector<double> X = jsonTrack["X"];
    std::vector<double> Y = jsonTrack["Y"];

    std::vector<double> X_inner = jsonTrack["X_i"];
    std::vector<double> Y_inner = jsonTrack["Y_i"];

    std::vector<double> X_outer = jsonTrack["X_o"];
    std::vector<double> Y_outer = jsonTrack["Y_o"];
//    TrackPos track_xy;
//    track_xy.X = Eigen::Map<Eigen::VectorXd>(X.data(), X.size());
//    track_xy.Y = Eigen::Map<Eigen::VectorXd>(Y.data(), Y.size());
    return {Eigen::Map<Eigen::VectorXd>(X.data(), X.size()), Eigen::Map<Eigen::VectorXd>(Y.data(), Y.size()),
            Eigen::Map<Eigen::VectorXd>(X_inner.data(), X_inner.size()), Eigen::Map<Eigen::VectorXd>(Y_inner.data(), Y_inner.size()),
            Eigen::Map<Eigen::VectorXd>(X_outer.data(), X_outer.size()), Eigen::Map<Eigen::VectorXd>(Y_outer.data(), Y_outer.size())};
}
}