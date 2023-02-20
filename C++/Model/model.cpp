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

#include "model.h"
namespace mpcc{
Model::Model()
:Ts_(1.0)
{
    std::cout << "default constructor, not everything is initialized properly" << std::endl;
}

Model::Model(double Ts,const PathToJson &path)
:Ts_(Ts),param_(Param(path.param_path))
{
    RK4_lib_ = (std::make_unique<CppAD::cg::LinuxDynamicLib<double>>(path.adcodegen_path+"/cppad_cg_RK4.so"));
    f_dyn_lib_ = (std::make_unique<CppAD::cg::LinuxDynamicLib<double>>(path.adcodegen_path+"/cppad_cg_f_dyn.so"));
    RK4_model_ = RK4_lib_->model("RK4");
    f_dyn_model_ = f_dyn_lib_->model("f_dyn");
}

StateVector Model::getF(const State &x,const Input &u) const
{
    return Eigen::Map<StateVector>((f_dyn_model_->ForwardZero(stateInputToVector(x,u))).data());
}

LinModelMatrix Model::discretizeModel(const State &x, const Input &u,const State &x_next) const
{
    // State x_lin = x;
    // x_lin.vxNonZero(param_.vx_zero);
    // std::vector<double> x_v_lin = stateInputToVector(x_lin,u);
    std::vector<double> x_v = stateInputToVector(x,u);
    Eigen::MatrixXd jac_vec = Eigen::Map<Eigen::Matrix<double,NX*(NX+NU),1>>((RK4_model_->Jacobian(x_v)).data());
    jac_vec.resize(NX+NU,NX);
    Eigen::MatrixXd jac = jac_vec.transpose();

    const A_MPC A_d = jac.block<NX,NX>(0,0);
    const B_MPC B_d = jac.block<NX,NU>(0,NX);

    StateVector x_RK4 = Eigen::Map<Eigen::Matrix<double,NX,1>>((RK4_model_->ForwardZero(x_v)).data());
    const g_MPC g_d =  -stateToVector(x_next) + x_RK4;
    return {A_d,B_d,g_d};
}

LinModelMatrix Model::getLinModel(const State &x, const Input &u, const State &x_next) const
{
    // discretize the system
    return discretizeModel(x,u,x_next);
}
}