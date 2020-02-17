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

#include "model_integrator_test.h"
namespace mpcc{
int testIntegrator(const PathToJson &path){
    double Ts = 0.02;
    const Integrator integrator = Integrator(Ts,path);

    // test integrator by comparing Euler forward to RK4
    // 3 differnet test points, hand picked, going straight and random

    //Integrator integrator;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Hand picked x and u
    StateVector error1;
    State xk1 = {0,0,0,2,0.1,-0.3,0.1,0.2,-0.1,2};
    Input uk1 = {0.2,-0.1,0};

    error1 = stateToVector(integrator.EF(xk1,uk1,Ts)) - stateToVector(integrator.RK4(xk1,uk1,Ts));
    std::cout << "hand picked point RK4 - EF error = " << error1.norm() << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // x and u corresponding to going straight with 1m/s
    StateVector error2;
    State xk2 = {0,0,0,1,0,0,0,0,0,0};
    Input uk2 ={0,0,0};

    error2 = stateToVector(integrator.EF(xk2,uk2,Ts)) - stateToVector(integrator.RK4(xk2,uk2,Ts));
    std::cout << "straight RK4 - EF error = " << error2.norm() << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // random x and u
    StateVector error3;
    StateVector xkr = StateVector::Random();
    InputVector ukr = InputVector::Random();
    State xk3 = vectorToState(xkr);
    Input uk3 = vectorToInput(ukr);

    error3 = stateToVector(integrator.EF(xk3,uk3,Ts)) - stateToVector(integrator.RK4(xk3,uk3,Ts));
    std::cout << "random RK4 - EF error = " <<  error3.norm() << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // test how good fit is
    if(error1.norm()/10.0 <= 0.3 && error2.norm()/10.0 <= 0.3 && error3.norm()/10.0 <= 0.3 ){
        return 1;
    }
    else{
        return 0;
    }

}

int testLinModel(const PathToJson &path){
    // test Liniear model by comparing it to RK4
    // 3 differnet test cases, hand picked, going straight and test how good linear model generalizes
    double Ts = 0.02;
    const Integrator integrator = Integrator(0.02,path);
    const Model model = Model(0.02,path);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Hand picked x and u
    StateVector error1;
    State xk1 = {0,0,0,2,0.1,-0.3,0.1,0.2,-0.1,2};
    Input uk1 = {0.2,-0.1,0};
    StateVector xk1_vec = stateToVector(xk1);
    InputVector uk1_vec = inputToVector(uk1);

    const LinModelMatrix lin_model_d1 = model.getLinModel(xk1,uk1);

    error1 = (lin_model_d1.A*xk1_vec + lin_model_d1.B*uk1_vec + lin_model_d1.g)  - stateToVector(integrator.RK4(xk1,uk1,Ts));
    std::cout << "hand picked point RK4 - lin error = " << error1.norm() << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // x and u corresponding to going straight with 1m/s
    StateVector error2;
    State xk2 = {0,0,0,1,0,0,0,0,0,0};
    Input uk2 = {0,0,0};
    StateVector xk2_vec = stateToVector(xk2);
    InputVector uk2_vec = inputToVector(uk2);

    const LinModelMatrix lin_model_d2 = model.getLinModel(xk2,uk2);

    error2 = (lin_model_d2.A*xk2_vec + lin_model_d2.B*uk2_vec + lin_model_d2.g)  - stateToVector(integrator.RK4(xk2,uk2,Ts));
    std::cout << "straight RK4 - lin error = " << error2.norm() << std::endl;
    //////////////////////////////////////////////////////////////////////////////////////////////////////////
    // generalization test
    // perturbe xk1 slightly, however us the model linearized around xk1 and uk1
    StateVector error3;
    State xk3;
    // xk3 is slightly perturbed version of xk1
    xk3 = xk1;
    xk3.vx += 0.2;  //vx
    xk3.vy += 0.05; //vy
    xk3.r += 0.8;  //r

    xk3.D += 0.2;  //D
    xk3.delta -= 0.05; //delta

    Input uk3;
    uk3 = uk1;
    //still linearize around xk1 and uk1
    StateVector xk3_vec = stateToVector(xk3);
    InputVector uk3_vec = inputToVector(uk3);
    error3 = (lin_model_d1.A*xk3_vec + lin_model_d1.B*uk3_vec + lin_model_d1.g)  - stateToVector(integrator.RK4(xk3,uk3,Ts));
//    std::cout << error3 << std::endl;
    std::cout << "generalization test RK4 - lin error = " << error3.norm() << std::endl;
    if(error1.norm()/10.0 <= 0.03 && error2.norm()/10.0 <= 0.03){
        return 1;
    }
    else{
        return 0;
    }
}
}