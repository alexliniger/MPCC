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

#include "arc_length_spline.h"

namespace mpcc{
ArcLengthSpline::ArcLengthSpline()
{ 
}
ArcLengthSpline::ArcLengthSpline(const PathToJson &path)
:param_(Param(path.param_path))
{
}

void ArcLengthSpline::setData(const Eigen::VectorXd &X_in,const Eigen::VectorXd &Y_in)
{
    // set input data if x and y have same length
    // compute arc length based on an piecewise linear approximation
    if(X_in.size() == Y_in.size()){
        path_data_.X = X_in;
        path_data_.Y = Y_in;
        path_data_.n_points = X_in.size();
        path_data_.s = compArcLength(X_in,Y_in);
    }
    else{
        std::cout << "input data does not have the same length" << std::endl;
    }
}

void ArcLengthSpline::setRegularData(const Eigen::VectorXd &X_in,const Eigen::VectorXd &Y_in,const Eigen::VectorXd &s_in) {
    // set final x-y data if x and y have same length
    // x-y points are space such that they are very close to arc length parametrized
    if(X_in.size() == Y_in.size()){
        path_data_.X = X_in;
        path_data_.Y = Y_in;
        path_data_.n_points = X_in.size();
        path_data_.s = s_in;
    }
    else{
        std::cout << "input data does not have the same length" << std::endl;
    }
}

Eigen::VectorXd ArcLengthSpline::compArcLength(const Eigen::VectorXd &X_in,const Eigen::VectorXd &Y_in) const
{
    //compute arc length based on straight line distance between the data points
    double dx,dy;
    double dist;

    int n_points = X_in.size();
    // initailize s as zero
    Eigen::VectorXd s;
    s.setZero(n_points);
//    std::cout << X_in << std::endl;
    for(int i=0;i<n_points-1;i++)
    {
        dx = X_in(i+1)-X_in(i);
        dy = Y_in(i+1)-Y_in(i);
        dist = std::sqrt(dx*dx + dy*dy);    //dist is straight line distance between points
        s(i+1) = s(i)+dist;       //s is cumulative sum of dist
    }
//    std::cout << s_in << std::endl;
    return s;
}

PathData ArcLengthSpline::resamplePath(const CubicSpline &initial_spline_x,const CubicSpline &initial_spline_y,const double total_arc_length) const
{
    // re-sample arc length parametrized X-Y spline path with N_spline data points
    // using equidistant arc length values
    // successively re-sample, computing the arc length and then fit the path should
    // result in close to equidistant points w.r.t. arc length

    // s -> "arc length" where points should be extracted
    // equilly spaced between 0 and current length of path
    PathData resampled_path;
    resampled_path.n_points=N_SPLINE;
    resampled_path.s.setLinSpaced(N_SPLINE,0,total_arc_length);

    // initialize new points as zero
    resampled_path.X.setZero(N_SPLINE);
    resampled_path.Y.setZero(N_SPLINE);

    // extract X-Y points
    for(int i=0;i<N_SPLINE;i++)
    {
        resampled_path.X(i) = initial_spline_x.getPoint(resampled_path.s(i));
        resampled_path.Y(i) = initial_spline_y.getPoint(resampled_path.s(i));
    }
    return resampled_path;
}

RawPath ArcLengthSpline::outlierRemoval(const Eigen::VectorXd &X_original,const Eigen::VectorXd &Y_original) const
{

    // remove points which are not at all equally spaced, to avoid fitting problems

    // compute mean distance between points and then process the points such that points
    // are not closer than 0.75 the mean distance

    double dx,dy;       // difference between points in x and y
    Eigen::VectorXd distVec;   // vector with all the distances
    double meanDist;    // mean distance
    double dist;        // temp variable for distance
    RawPath resampled_path;
    int k = 0;          // indecies
    int j = 0;

    if (X_original.size() != Y_original.size()){
        //error
    }
//    std::cout << X_original << std::endl;

    int n_points = X_original.size();

    // initialize with zero
    resampled_path.X.setZero(n_points);
    resampled_path.Y.setZero(n_points);



    // compute distance between points in X-Y data
    distVec.setZero(n_points-1);
    for(int i=0;i<n_points-1;i++){
        dx = X_original(i+1)-X_original(i);
        dy = Y_original(i+1)-Y_original(i);
        distVec(i) = std::sqrt(dx*dx + dy*dy);
    }
    // compute mean distance between points
    meanDist = distVec.sum()/(n_points-1);

    // compute the new points
    // start point is the original start point
    resampled_path.X(k) = X_original(k);
    resampled_path.Y(k) = Y_original(k);
    k++;
    for(int i=1;i<n_points-1;i++){
        // compute distance between currently checked point and the one last added to the new X-Y path
        dx = X_original(i)-X_original(j);
        dy = Y_original(i)-Y_original(j);
        dist = std::sqrt(dx*dx + dy*dy);
        // if this distance is smaller than 0.7 the mean distance add this point to the new X-Y path
        if(dist >= 0.7*meanDist)
        {
            resampled_path.X(k) = X_original(i);
            resampled_path.Y(k) = Y_original(i);
            k++;
            j = i;
        }
    }
    // always add the last point
    resampled_path.X(k) = X_original(n_points-1);
    resampled_path.Y(k) = Y_original(n_points-1);
    k++;

//    std::cout << "not resiszed " << X_new.transpose() << std::endl;
    // set the new X-Y data
//    setData(X.head(k),Y.head(k));
    resampled_path.X.conservativeResize(k);
    resampled_path.Y.conservativeResize(k);

//    std::cout << "resiszed " << X_new.transpose()  << std::endl;
    return resampled_path;
}

double ArcLengthSpline::unwrapInput(double x) const
{
    double x_max = getLength();
    return x - x_max*std::floor(x/x_max);
}

void ArcLengthSpline::fitSpline(const Eigen::VectorXd &X,const Eigen::VectorXd &Y)
{
    // successively fit spline -> re-sample path -> compute arc length
    // temporary spline class only used for fitting
    Eigen::VectorXd s_approximation;
    PathData first_refined_path,second_refined_path;
    double total_arc_length;

    s_approximation = compArcLength(X,Y);
//    std::cout << s_approximation << std::endl;
    total_arc_length = s_approximation(s_approximation.size()-1);

    CubicSpline first_spline_x,first_spline_y;
    CubicSpline second_spline_x,second_spline_y;
    // 1. spline fit
    first_spline_x.genSpline(s_approximation,X,false);
    first_spline_y.genSpline(s_approximation,Y,false);
    // 1. re-sample
    first_refined_path = resamplePath(first_spline_x,first_spline_y,total_arc_length);
    s_approximation = compArcLength(first_refined_path.X,first_refined_path.Y);

    total_arc_length = s_approximation(s_approximation.size()-1);
    ////////////////////////////////////////////
    // 2. spline fit
    second_spline_x.genSpline(s_approximation,first_refined_path.X,false);
    second_spline_y.genSpline(s_approximation,first_refined_path.Y,false);
    // 2. re-sample
    second_refined_path = resamplePath(second_spline_x,second_spline_y,total_arc_length);
    ////////////////////////////////////////////
    setRegularData(second_refined_path.X,second_refined_path.Y,second_refined_path.s);
//    setData(second_refined_path.X,second_refined_path.Y);
    // Final spline fit with fixed Delta_s
    spline_x_.genSpline(path_data_.s,path_data_.X,true);
    spline_y_.genSpline(path_data_.s,path_data_.Y,true);


}

void ArcLengthSpline::gen2DSpline(const Eigen::VectorXd &X,const Eigen::VectorXd &Y)
{
    // generate 2-D arc length parametrized spline given X-Y data

    // remove outliers, depending on how iregular the points are this can help
    RawPath clean_path;
    clean_path = outlierRemoval(X,Y);
    // successively fit spline and re-sample
    fitSpline(clean_path.X,clean_path.Y);

}


Eigen::Vector2d ArcLengthSpline::getPostion(const double s) const
{
    Eigen::Vector2d s_path;
    s_path(0) = spline_x_.getPoint(s);
    s_path(1) = spline_y_.getPoint(s);

    return s_path;
}

Eigen::Vector2d ArcLengthSpline::getDerivative(const double s) const
{
    Eigen::Vector2d ds_path;
    ds_path(0) = spline_x_.getDerivative(s);
    ds_path(1) = spline_y_.getDerivative(s);

    return ds_path;
}

Eigen::Vector2d ArcLengthSpline::getSecondDerivative(const double s) const
{
    Eigen::Vector2d dds_path;
    dds_path(0) = spline_x_.getSecondDerivative(s);
    dds_path(1) = spline_y_.getSecondDerivative(s);

    return dds_path;
}

double ArcLengthSpline::getLength() const
{
    return path_data_.s(path_data_.n_points-1);
}

double ArcLengthSpline::porjectOnSpline(const State &x) const
{
    Eigen::Vector2d pos;
    pos(0) = x.X;
    pos(1) = x.Y;
    double s_guess = x.s;
    Eigen::Vector2d pos_path = getPostion(s_guess);

    double s_opt = s_guess;
    double dist = (pos-pos_path).norm();

    if (dist >= param_.max_dist_proj)
    {
        std::cout << "dist too large" << std::endl;
        Eigen::ArrayXd diff_x_all = path_data_.X.array() - pos(0);
        Eigen::ArrayXd diff_y_all = path_data_.Y.array() - pos(1);
        Eigen::ArrayXd dist_square = diff_x_all.square() + diff_y_all.square();
        std::vector<double> dist_square_vec(dist_square.data(),dist_square.data() + dist_square.size());
        auto min_iter = std::min_element(dist_square_vec.begin(),dist_square_vec.end());
        s_opt = path_data_.s(std::distance(dist_square_vec.begin(), min_iter));
    }
    double s_old = s_opt;
    for(int i=0; i<20; i++)
    {
        pos_path = getPostion(s_opt);
        Eigen::Vector2d ds_path = getDerivative(s_opt);
        Eigen::Vector2d dds_path = getSecondDerivative(s_opt);
        Eigen::Vector2d diff = pos_path - pos;
        double jac = 2.0 * diff(0) * ds_path(0) + 2.0 * diff(1) * ds_path(1);
        double hessian = 2.0 * ds_path(0) * ds_path(0) + 2.0 * diff(0) * dds_path(0) +
                         2.0 * ds_path(1) * ds_path(1) + 2.0 * diff(1) * dds_path(1);
        // Newton method
        s_opt -= jac/hessian;
        s_opt = unwrapInput(s_opt);

//        std::cout << std::abs(s_old - s_opt) << std::endl;
        if(std::abs(s_old - s_opt) <= 1e-5)
            return s_opt;
        s_old = s_opt;
    }
    // something is strange if it did not converge within 20 iterations, give back the initial guess
    return s_guess;
}
}