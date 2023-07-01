//
// Created by alexliniger on 03.11.20.
//

#include "boost_splines.h"

namespace mpcc{
BoostSplines::BoostSplines()
{
}

BoostSplines::BoostSplines(const PathToJson &path)
: param_(Param(path.param_path))
{
}

void BoostSplines::genSplines(const TrackFull &track)
{
    genPathSplines(track.s,track.X,track.Y);
    genBoarderSplines(track.s,track.n_left,track.n_right);
    genVeloSplines(track.s,track.v);
}

void BoostSplines::genPathSplines(const Eigen::VectorXd &s,const Eigen::VectorXd &X, const Eigen::VectorXd &Y)
{
    double diff_s = s[1] - s[0];
    boost::math::interpolators::cardinal_cubic_b_spline<double> temp_x(X.begin(), X.end(), 0.0, diff_s);
    boost::math::interpolators::cardinal_cubic_b_spline<double> temp_y(Y.begin(), Y.end(), 0.0, diff_s);

    spline_x_ = temp_x;
    spline_y_ = temp_y;

    path_data_.X = Eigen::Map<const Eigen::VectorXd>(X.data(), X.size());
    path_data_.Y = Eigen::Map<const Eigen::VectorXd>(Y.data(), Y.size());
    path_data_.s = Eigen::Map<const Eigen::VectorXd>(s.data(), s.size());
    path_data_.n_points = int(s.size());
    length_s_ = s.tail(1)(0);
}

void BoostSplines::genBoarderSplines(const Eigen::VectorXd &s,const Eigen::VectorXd &n_left, const Eigen::VectorXd &n_right)
{
    double diff_s = s[1] - s[0];
    boost::math::interpolators::cardinal_cubic_b_spline<double> temp_n_left(n_left.begin(), n_left.end(), 0.0, diff_s);
    boost::math::interpolators::cardinal_cubic_b_spline<double> temp_n_right(n_right.begin(), n_right.end(), 0.0, diff_s);

    spline_n_left_ = temp_n_left;
    spline_n_right_ = temp_n_right;
}

void BoostSplines::genVeloSplines(const Eigen::VectorXd &s,const Eigen::VectorXd &v)
{
    double diff_s = s[1] - s[0];
    boost::math::interpolators::cardinal_cubic_b_spline<double> temp_v(v.begin(), v.end(), 0.0, diff_s);

    spline_v_ = temp_v;
}

Eigen::Vector2d BoostSplines::getPostion(double s) const
{
    double X = spline_x_(unwrapInput(s));
    double Y = spline_y_(unwrapInput(s));
    return {X,Y};
}

Eigen::Vector2d BoostSplines::getDerivative(const double s) const
{
    double ds_X = spline_x_.prime(unwrapInput(s));
    double ds_Y = spline_y_.prime(unwrapInput(s));
    return {ds_X,ds_Y};
}

Eigen::Vector2d BoostSplines::getSecondDerivative(const double s) const
{
    double dds_X = spline_x_.double_prime(unwrapInput(s));
    double dds_Y = spline_y_.double_prime(unwrapInput(s));
    return {dds_X,dds_Y};
}

double BoostSplines::getNLeft(double s) const
{
    return spline_n_left_(unwrapInput(s));
}

double BoostSplines::getNRight(double s) const
{
    return spline_n_right_(unwrapInput(s));
}

double BoostSplines::getVelocity(double s) const
{
    return spline_v_(unwrapInput(s));
}

double BoostSplines::getLength() const
{
    return length_s_;
}

double BoostSplines::projectOnSpline(const State &x) const
{
    Eigen::Vector2d pos;
    pos(0) = x.X;
    pos(1) = x.Y;
    double s_guess = x.s;
    Eigen::Vector2d pos_path = {spline_x_(s_guess),spline_y_(s_guess)};

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
        Eigen::Vector2d ds_path = {spline_x_.prime(s_opt),spline_y_.prime(s_opt)};
        Eigen::Vector2d dds_path = {spline_x_.double_prime(s_opt),spline_y_.double_prime(s_opt)};
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

double BoostSplines::unwrapInput(double x) const{
    double x_max = getLength();
    return x - x_max*std::floor(x/x_max);
}

}