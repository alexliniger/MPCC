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

#ifndef MPCC_TRACK_H
#define MPCC_TRACK_H

#include "config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <nlohmann/json.hpp>

namespace mpcc {
//used namespace
using json = nlohmann::json;

struct TrackPos {
    const Eigen::VectorXd X;
    const Eigen::VectorXd Y;

    const Eigen::VectorXd X_inner;
    const Eigen::VectorXd Y_inner;

    const Eigen::VectorXd X_outer;
    const Eigen::VectorXd Y_outer;
};

 struct TrackFull {
    const Eigen::VectorXd X;
    const Eigen::VectorXd Y;

    const Eigen::VectorXd X_inner;
    const Eigen::VectorXd Y_inner;

    const Eigen::VectorXd X_outer;
    const Eigen::VectorXd Y_outer;

    const Eigen::VectorXd n_left;
    const Eigen::VectorXd n_right;

    const Eigen::VectorXd s;
    const Eigen::VectorXd v;
 };


class Track {
public:
    Track(std::string file);
    TrackPos getTrack() const;
    TrackFull getTrackFull() const;

private:
    Eigen::VectorXd X_;
    Eigen::VectorXd Y_;

    Eigen::VectorXd X_inner_;
    Eigen::VectorXd Y_inner_;

    Eigen::VectorXd X_outer_;
    Eigen::VectorXd Y_outer_;

    Eigen::VectorXd n_left_;
    Eigen::VectorXd n_right_;

    Eigen::VectorXd s_;

    Eigen::VectorXd v_;
};
};

#endif //MPCC_TRACK_H
