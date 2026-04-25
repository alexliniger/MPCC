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

#include "spline_test.h"
namespace mpcc{
int testSpline() {

    //Fit spline to a cos and test if the spine interpolation is accurate
    //Also test if first and second derivatives are accurate

    CubicSpline spline;

    Eigen::VectorXd x;
    Eigen::VectorXd y;

    int NT = 50;    //number of "training" points
    int NV = 100;   //number of validation points
    x.setLinSpaced(NT,0,2*M_PI);
    y.setZero(NT,1);
//    std::cout << x.size() << std::endl;
    //set y ot cos(x)
    for(int i=0;i<x.size();i++){
        y(i) = std::cos(x(i));
    }
    //give data to spline class and compute spline parameters
    spline.genSpline(x,y,true);

    //spline test outputs and true outputs -> initialize to zero
    Eigen::VectorXd xt;
    Eigen::VectorXd yt,ytt;
    Eigen::VectorXd dyt,dytt;
    Eigen::VectorXd ddyt,ddytt;

    xt.setLinSpaced(NV,0,2*M_PI); //test x points
    yt.setZero(NV,1);
    ytt.setZero(NV,1);

    dyt.setZero(NV,1);
    dytt.setZero(NV,1);

    ddyt.setZero(NV,1);
    ddytt.setZero(NV,1);

    //compute spline outputs and true outputs given the y = cos(x)
    for(int i=0;i<xt.size();i++){
        yt(i) = spline.getPoint(xt(i));
        ytt(i) = std::cos(xt(i));

        dyt(i) = spline.getDerivative(xt(i));
        dytt(i) = -std::sin(xt(i));

        ddyt(i) = spline.getSecondDerivative(xt(i));
        ddytt(i) = -std::cos(xt(i));
    }

    // output how good the fit is
    std::cout << "mean fit error: " << (yt-ytt).sum()/NV << std::endl;

    std::cout << "mean fit derivative error: " << (dyt-dytt).sum()/NV << std::endl;

    std::cout << "mean fit second derivative error: " << (ddyt-ddytt).sum()/NV << std::endl;

    //return true if fit is good enough
    if ((yt-ytt).sum()/NV <= 1e-4 && (dyt-dytt).sum()/NV <= 1e-3 && (ddyt-ddytt).sum()/NV <= 1e-1){
        return 1;
    }
    else {
        return 0;
    }

}

int testArcLengthSpline(const PathToJson &path){
    // test 2-D arc length spline approach
    // given a circle with randomly distributed points

    ArcLengthSpline twoDspline = ArcLengthSpline(path);

    int NT = 50;    //number of "training" points
    int NV = 200;   //number of validation points points

    Eigen::VectorXd X;
    Eigen::VectorXd Y;
    Eigen::VectorXd phiR;

    X.setZero(NT);
    Y.setZero(NT);
    // randomly distribute training points around circle
    // generate random agles between [0,2pi]
    phiR.setRandom(NT);
    phiR = (phiR/2.0+0.5*Eigen::VectorXd::Ones(NT))*2*M_PI;
    // sort training points
    std::sort(phiR.data(), phiR.data()+phiR.size());
    // fix start and end point
    phiR(0) = 0;
    phiR(NT-1) = 2*M_PI;
    // compute circle points
    for(int i=0;i<NT;i++){
        X(i) = std::cos(phiR(i));
        Y(i) = std::sin(phiR(i));
    }

    // give points to arc length based 2-D curve fitting
    twoDspline.gen2DSpline(X,Y);

    //validation points
    Eigen::VectorXd Xv;
    Eigen::VectorXd Yv;
    Eigen::VectorXd phiv;
    //error between true circle and spline approximation
    Eigen::VectorXd error;
    // uniformly space validaiton points between [0,2pi]
    Xv.setZero(NV);
    Yv.setZero(NV);
    phiv.setLinSpaced(NV,0,2*M_PI);

    error.setZero(NV);

    Eigen::Vector2d pos;

    for(int i=0;i<NV;i++){
        pos = twoDspline.getPostion(phiv(i));
        error(i) = std::sqrt(std::pow(pos(0) - std::cos(phiv(i)),2) + std::pow(pos(1) -std::sin(phiv(i)),2));
    }
    std::cout << "norm of error = " << error.norm() << std::endl;

    if (error.norm()/(double)NV <= 0.03)
        return 1;
    else
        return 0;

    //    ArcLengthSpline two_d_spline;
//    two_d_spline.gen2DSpline(track_xy.X,track_xy.Y);
//    Eigen::VectorXd s;
//    s.setLinSpaced(1000,0,two_d_spline.getLength());
//    Eigen::Vector2d pos;
//    std::vector<double> plot_x;
//    std::vector<double> plot_y;
//    std::vector<double> kappa;
//    for(int i=0;i<1000;i++){
//        std::cout << s(i) << std::endl;
//        pos = two_d_spline.getPostion(s(i));
//        plot_x.push_back(pos(0));
//        plot_y.push_back(pos(1));
//        // reference path derivatives
//        const Eigen::Vector2d dpos_ref = two_d_spline.getDerivative(s(i));
//        const double dx_ref = dpos_ref(0);
//        const double dy_ref = dpos_ref(1);
//        // angle of the reference path
//        const double theta_ref = atan2(dy_ref,dx_ref);
//        // second order derivatives
//        Eigen::Vector2d ddpos_ref = two_d_spline.getSecondDerivative(s(i));
//        const double ddx_ref = ddpos_ref(0);
//        const double ddy_ref = ddpos_ref(1);
//        // curvature
//        double dtheta_ref = 0.0;
//        if(std::abs(ddx_ref*ddx_ref + ddy_ref*ddy_ref)>= 1e0){
//            dtheta_ref = (dx_ref*ddy_ref - dy_ref*ddx_ref)
//                         /(ddx_ref*ddx_ref + ddy_ref*ddy_ref); //curvature
//        }
//        kappa.push_back(dtheta_ref);
//    }
//    plt::plot(plot_x,plot_y);
//    plt::show();
//
//    plt::plot(kappa);
//    plt::show();
}
}