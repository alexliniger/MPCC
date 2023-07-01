//
// Created by alexliniger on 03.11.20.
//

#ifndef MPCC_BOOST_SPLINES_H
#define MPCC_BOOST_SPLINES_H

#include "types.h"
#include "config.h"
#include "Params/params.h"
#include "Params/track.h"
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
// data struct
struct PathData{
    Eigen::VectorXd X;
    Eigen::VectorXd Y;
    Eigen::VectorXd s;
    int n_points;
};

namespace mpcc{
class BoostSplines{
public:
    void genSplines(const TrackFull &track);

    Eigen::Vector2d getPostion(double s) const;
    Eigen::Vector2d getDerivative(double s) const;
    Eigen::Vector2d getSecondDerivative(double s) const;
    double getNLeft(double s) const;
    double getNRight(double s) const;
    double getVelocity(double s) const;
    double getLength() const;
    double projectOnSpline(const State &x) const;

    BoostSplines();
    BoostSplines(const PathToJson &path);
private:
    void genPathSplines(const Eigen::VectorXd &s,const Eigen::VectorXd &X, const Eigen::VectorXd &Y);
    void genBoarderSplines(const Eigen::VectorXd &s,const Eigen::VectorXd &n_left, const Eigen::VectorXd &n_right);
    void genVeloSplines(const Eigen::VectorXd &s,const Eigen::VectorXd &v);
    double unwrapInput(double x) const;

    double length_s_;
    PathData path_data_;      // initial data and data used for successive fitting
    boost::math::interpolators::cardinal_cubic_b_spline<double> spline_x_;
    boost::math::interpolators::cardinal_cubic_b_spline<double> spline_y_;
    boost::math::interpolators::cardinal_cubic_b_spline<double> spline_n_left_;
    boost::math::interpolators::cardinal_cubic_b_spline<double> spline_n_right_;
    boost::math::interpolators::cardinal_cubic_b_spline<double> spline_v_;

    Param param_;
};
}

#endif //MPCC_BOOST_SPLINES_H
