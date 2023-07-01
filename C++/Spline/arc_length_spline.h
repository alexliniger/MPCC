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

#ifndef MPCC_ARC_LENGTH_SPLINE_H
#define MPCC_ARC_LENGTH_SPLINE_H

#include "cubic_spline.h"
#include "types.h"
#include "Params/params.h"
#include <map>

namespace mpcc{
//return value
struct RawPath{
    Eigen::VectorXd X;
    Eigen::VectorXd Y;
};
// data struct
struct PathData{
    Eigen::VectorXd X;
    Eigen::VectorXd Y;
    Eigen::VectorXd s;
    int n_points;
};

class ArcLengthSpline {
public:
    // X and Y spline used for final spline fit
    void gen2DSpline(const Eigen::VectorXd &X,const Eigen::VectorXd &Y);
    Eigen::Vector2d getPostion(double) const;
    Eigen::Vector2d getDerivative(double) const;
    Eigen::Vector2d getSecondDerivative(double) const;
    double getLength() const;
    double projectOnSpline(const State &x) const;

    ArcLengthSpline();
    ArcLengthSpline(const PathToJson &path);
    // void setParam(const Param &param) { param_ = param; };

private:
    void setData(const Eigen::VectorXd &X_in,const Eigen::VectorXd &Y_in);
    void setRegularData(const Eigen::VectorXd &X_in,const Eigen::VectorXd &Y_in,const Eigen::VectorXd &s_in);
    Eigen::VectorXd compArcLength(const Eigen::VectorXd &X_in,const Eigen::VectorXd &Y_in) const;
    PathData resamplePath(const CubicSpline &initial_spline_x,const CubicSpline &initial_spline_y,double total_arc_length) const;
    RawPath outlierRemoval(const Eigen::VectorXd &X_original,const Eigen::VectorXd &Y_original) const;
    void fitSpline(const Eigen::VectorXd &X,const Eigen::VectorXd &Y);
    double unwrapInput(double x) const;

    PathData path_data_;      // initial data and data used for successive fitting
//    PathData pathDataFinal; // final data
    CubicSpline spline_x_;
    CubicSpline spline_y_;
    Param param_;
};
}
#endif //MPCC_ARC_LENGTH_SPLINE_H