#include <iosfwd>
#include <iostream>
#include <memory>

#include "ad_dynamics.h"
#include <iosfwd>
#include <vector>
#include <cppad/cg.hpp>
#include "config.h"
#include "types.h"
#include "Params/params.h"
using json = nlohmann::json;

int main(void) {
    /***************************************************************************
     *                       Use the dynamic library
     **************************************************************************/
    using namespace mpcc;
    std::ifstream iConfig("../../Params/config.json");
    json jsonConfig;
    iConfig >> jsonConfig;

    ADDynamics ad_dyn(jsonConfig["Ts"],std::string("../../Params/model.json"));
    ad_dyn.genLibraryRK4();
    ad_dyn.genLibraryGetF();
    ad_dyn.genLibraryTireFront();
    ad_dyn.genLibraryTireRear();

    std::unique_ptr<CppAD::cg::LinuxDynamicLib<double>> myLib = std::make_unique<CppAD::cg::LinuxDynamicLib<double>>("cppad_cg_RK4.so");
    std::unique_ptr<CppAD::cg::GenericModel<double>> model;
    model = myLib->model("RK4");

    Eigen::Vector<double,15> x;
    x.setZero();
    x(3) = 10.0;
    std::vector<double> xv(x.data(), x.data() + x.size());
    Eigen::MatrixXd jac = Eigen::Map<Eigen::Matrix<double,NX*(NX+NU),1>>((model->Jacobian(xv)).data());
    jac.resize(NX+NU,NX);
    Eigen::MatrixXd jac2 = jac.transpose();
    std::cout << jac2 << std::endl;
    StateVector fw = Eigen::Map<StateVector>((model->ForwardZero(xv)).data());
    std::cout << std::endl;
    std::cout << fw << std::endl;
    std::cout << std::endl;


    std::unique_ptr<CppAD::cg::LinuxDynamicLib<double>> cont_dyn_lib = std::make_unique<CppAD::cg::LinuxDynamicLib<double>>("cppad_cg_f_dyn.so");
    std::unique_ptr<CppAD::cg::GenericModel<double>> model_cont_dyn;
    model_cont_dyn = cont_dyn_lib->model("f_dyn");
    StateVector f = Eigen::Map<StateVector>((model_cont_dyn->ForwardZero(xv)).data());
    std::cout << f << std::endl;

    Eigen::Vector<double,NX> s;
    s.setZero();
    s(3) = 10.0;
    s(5) = 1.0;
    std::vector<double> sv(s.data(), s.data() + s.size());

    std::unique_ptr<CppAD::cg::LinuxDynamicLib<double>> tcf_lib = std::make_unique<CppAD::cg::LinuxDynamicLib<double>>("cppad_cg_TCR.so");
    std::unique_ptr<CppAD::cg::GenericModel<double>> model_tcf;
    model_tcf = tcf_lib->model("TireConRear");
    std::vector<double> c_0 = model_tcf->ForwardZero(sv);
    std::cout << c_0[0] << std::endl;
    C_i_MPC C_tire_front_constraint = Eigen::Map<C_i_MPC>((model_tcf->Jacobian(sv)).data());
    std::cout << C_tire_front_constraint << std::endl;
}