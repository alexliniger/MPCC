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

#include "constratins_test.h"
namespace mpcc{
void genRoundTrack(ArcLengthSpline &track){
    int NT = 100;    //number of track points
    double TrackRadius = 1.0; // track radius

    Eigen::VectorXd X;
    Eigen::VectorXd Y;
    Eigen::VectorXd phiR;

    X.setZero(NT);
    Y.setZero(NT);
//    // randomly distribute training points around circle
//    // generate random agles between [0,2pi]
//    phiR.setRandom(NT);
//    phiR = (phiR/2.0+0.5*Eigen::VectorXd::Ones(NT))*2*M_PI;
//    // sort training points
//    std::sort(phiR.data(), phiR.data()+phiR.size());
//    // fix start and end point
//    phiR(0) = 0;
//    phiR(NT-1) = 2*M_PI;
    // generate equally spaced points around circle
    phiR.setLinSpaced(NT,0,2*M_PI);
    // compute circle points
    for(int i=0;i<NT;i++){
        X(i) = TrackRadius*std::cos(phiR(i));
        Y(i) = TrackRadius*std::sin(phiR(i));
    }

    // give points to arc length based 2-D curve fitting
    track.gen2DSpline(X,Y);



//    std::vector<double> Xv(NT),Yv(NT),Xsv(NT*10),Ysv(NT*10);
//    for(int i=0;i<NT;i++){
//        Xv[i] = X(i);
//        Yv[i] = Y(i);
//    }
//    Eigen::VectorXd s;
//    s.setLinSpaced(NT*10,0,2*M_PI);
//    for(int i=0;i<10*NT;i++){
//        Xsv[i] = track.splineX.ppval(s(i));
//        Ysv[i] = track.splineY.ppval(s(i));
//    }
//
//    plt::plot(Xv,Yv,"--g");
//    plt::draw();
//    plt::plot(Xsv,Ysv,"--r");
//    plt::draw();
//    plt::show();

}

int testAlphaConstraint(const PathToJson &path){
    Constraints constraints = Constraints(0.02,path);
    ArcLengthSpline track = ArcLengthSpline(path);
    Model model = Model(0.02,path);
    Param param = Param(path.param_path);
    // track.setParam(param);
    
    double alpha_f;
    genRoundTrack(track);

    StateVector xk0_vec, xk1_vec, xk2_vec;
    xk0_vec << 0,0,0,2,0,0,0,0,0,0;
    xk1_vec << 0,0,0,2,-0.8,-3,0.0,0.2,0.1,2;
    xk2_vec  << 0,0,0,2,-0.8,-3,0.0,0.2,0.01,2;
    InputVector uk0_vec, uk1_vec, uk2_vec;
    uk0_vec << 0,0,0;
    uk1_vec << 0.2,-0.1,0;
    uk1_vec << 0.2,-0.1,0;

    State xk0 = vectorToState(xk0_vec);
    Input uk0 = vectorToInput(uk0_vec);
    State xk1 = vectorToState(xk1_vec);
    Input uk1 = vectorToInput(uk1_vec);
    State xk2 = vectorToState(xk2_vec);
    Input uk2 = vectorToInput(uk2_vec);


    alpha_f = model.getSlipAngleFront(xk1);

    ConstrainsMatrix constraints_0 = constraints.getConstraints(track,xk0,uk0);
    std::cout << "true alpha " << alpha_f  << std::endl;
    std::cout << constraints_0.dl(2) << "<=" << constraints_0.C.row(2)*xk1_vec <<  "<=" << constraints_0.du(2) << std::endl;

    if ((constraints_0.dl(2)<= constraints_0.C.row(2)*xk1_vec && constraints_0.C.row(2)*xk1_vec  <= constraints_0.du(2)) ^ (-param.max_alpha <= alpha_f && alpha_f  <= param.max_alpha)){
        return  0;
    }

    ConstrainsMatrix constraints_1 = constraints.getConstraints(track,xk1,uk1);
    std::cout << constraints_1.dl(2) << "<=" << constraints_1.C.row(2)*xk1_vec <<  "<=" << constraints_1.du(2) << std::endl;

    if ((constraints_1.dl(2)<= constraints_1.C.row(2)*xk1_vec && constraints_1.C.row(2)*xk1_vec  <= constraints_1.du(2)) ^ (-param.max_alpha <= alpha_f && alpha_f  <= param.max_alpha)){
        return  0;
    }


    alpha_f = model.getSlipAngleFront(xk2);

    ConstrainsMatrix constraints_2 = constraints.getConstraints(track,xk2,uk2);

    std::cout << "true alpha " << alpha_f  << std::endl;
    std::cout << constraints_2.dl(2) << "<=" << constraints_2.C.row(2)*xk2_vec <<  "<=" << constraints_2.du(2) << std::endl;

    if ((constraints_2.dl(2)<= constraints_2.C.row(2)*xk2_vec && constraints_2.C.row(2)*xk2_vec  <= constraints_2.du(2)) ^ (-param.max_alpha <= alpha_f && alpha_f  <= param.max_alpha)){
        return  0;
    }

    return 1;

}


int testTireForceConstraint(const PathToJson &path) {
    Constraints constraints = Constraints(0.02,path);
    ArcLengthSpline track = ArcLengthSpline(path);
    Model model = Model(0.02,path);
    Param param = Param(path.param_path);


    double Frx, Fry;
    double tireForce;
    double maxForce;

    genRoundTrack(track);

    StateVector xk0_vec, xk1_vec, xk2_vec;
    xk0_vec << 0, 0, 0, 2, 0, 0, 0, 0, 0, 0;
    xk1_vec << 0, 0, 0, 1, -0.6, -3.5, 0.0, 0.6, 0.1, 2;
    xk2_vec << 0, 0, 0, 1, -0.8, -3, 0.0, 0.5, 0.01, 2;
    InputVector uk0_vec, uk1_vec, uk2_vec;
    uk0_vec << 0, 0, 0;
    uk1_vec << 0.2, -0.1, 0;
    uk2_vec << 0, 0, 0;

    State xk0 = vectorToState(xk0_vec);
    Input uk0 = vectorToInput(uk0_vec);
    State xk1 = vectorToState(xk1_vec);
    Input uk1 = vectorToInput(uk1_vec);
    State xk2 = vectorToState(xk2_vec);
    Input uk2 = vectorToInput(uk2_vec);

    TireForces f_rear = model.getForceRear(xk1);

    tireForce = std::sqrt(std::pow(param.e_long*f_rear.F_x,2) + std::pow(f_rear.F_y,2));
    maxForce = param.e_eps*param.Dr;

    ConstrainsMatrix constraints_mat = constraints.getConstraints(track,xk2,uk2);

    //    std::cout << C.row(1) <<   std::endl;
    std::cout << 0 << "<=" << tireForce <<  "<=" << maxForce << std::endl;
    std::cout << constraints_mat.dl(1) << "<=" << constraints_mat.C.row(1)*xk1_vec <<  "<=" << constraints_mat.du(1) << std::endl;

    return  1;
}


int testTrackConstraint(const PathToJson &path) {
    Constraints constraints = Constraints(0.02,path);
    ArcLengthSpline track = ArcLengthSpline(path);

    genRoundTrack(track);

    StateVector xk0_vec, xk1_vec, xk2_vec;
    xk0_vec << 1, 0, 0, 2, 0, 0, 0, 0, 0, 0;
    xk1_vec << 1.2, 0, 0, 2, 0, 0, 0, 0, 0, 0;
    xk2_vec << 0.8, 0, 0, 2, 0, 0, 0, 0, 0, 0;
    InputVector uk0_vec;
    uk0_vec << 0, 0, 0;

    State xk0 = vectorToState(xk0_vec);
    Input uk0 = vectorToInput(uk0_vec);
    State xk1 = vectorToState(xk1_vec);
    State xk2 = vectorToState(xk2_vec);

    ConstrainsMatrix constraints_mat = constraints.getConstraints(track,xk0,uk0);

//    std::cout << "inside track case " << dl(0) << "<=" << C.row(0)*xk0 <<  "<=" << du(0) << std::endl;
//    std::cout << "right side case "   << dl(0) << "<=" << C.row(0)*xk1 <<  "<=" << du(0) << std::endl;
//    std::cout << "left side case "    << dl(0) << "<=" << C.row(0)*xk2 <<  "<=" << du(0) << std::endl;

    if (!(constraints_mat.dl(0)<= constraints_mat.C.row(0)*xk0_vec && constraints_mat.C.row(0)*xk0_vec  <= constraints_mat.du(0))){
        return  0;
    }
    if (constraints_mat.dl(0)<= constraints_mat.C.row(0)*xk1_vec && constraints_mat.C.row(0)*xk1_vec  <= constraints_mat.du(0)){
        return  0;
    }
    if (constraints_mat.dl(0)<= constraints_mat.C.row(0)*xk2_vec && constraints_mat.C.row(0)*xk2_vec  <= constraints_mat.du(0)){
        return  0;
    }

    return  1;
}
}