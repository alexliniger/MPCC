//
// Created by alexliniger on 24.07.20.
//

#include "ad_dynamics.h"
namespace mpcc{
ADDynamics::ADDynamics():Ts_(1.0)
{
    std::cout << "default constructor, not everything is initialized properly" << std::endl;
}
ADDynamics::ADDynamics(double Ts,const std::string &path)
:Ts_(Ts),param_(path)
{
}
void ADDynamics::genLibraryRK4() {
    // independent variable vector
    std::vector<ADCG> z(NX + NU);
    for(int i = 0; i< NX+NU;i++)
        z[i] = 0.0;
    z[3] = 10.;
    Independent(z);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // RK4
    // dependent variable vector
    std::vector<ADCG> x_plus(NX);
    // the model equation
    x_plus = RK4(z);
    ADFun<CGD> fun(z, x_plus);

    // generates source code
    ModelCSourceGen<double> cgen(fun, "RK4");
    cgen.setCreateJacobian(true);
    ModelLibraryCSourceGen<double> libcgen(cgen);

    // compile source code
    DynamicModelLibraryProcessor<double> p(libcgen, "cppad_cg_RK4");

    GccCompiler<double> compiler;
    std::unique_ptr<DynamicLib<double>> dynamicLibRK4 = p.createDynamicLibrary(compiler);

    // save to files (not really required)
    SaveFilesModelLibraryProcessor<double> p2(libcgen);
    p2.saveSources();
}

void ADDynamics::genLibraryGetF() {
    // independent variable vector
    std::vector<ADCG> z(NX + NU);
    for(int i = 0; i< NX+NU;i++)
        z[i] = 0.0;
    z[3] = 10.;
    Independent(z);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // RK4
    // dependent variable vector
    std::vector<ADCG> x_dot(NX);
    // the model equation
    x_dot = f_dyn(z);
    ADFun<CGD> fun(z, x_dot);

    // generates source code
    ModelCSourceGen<double> cgen(fun, "f_dyn");
    cgen.setCreateJacobian(true);
    ModelLibraryCSourceGen<double> libcgen(cgen);

    // compile source code
    DynamicModelLibraryProcessor<double> p(libcgen, "cppad_cg_f_dyn");

    GccCompiler<double> compiler;
    std::unique_ptr<DynamicLib<double>> dynamicLibRK4 = p.createDynamicLibrary(compiler);

    // save to files (not really required)
    SaveFilesModelLibraryProcessor<double> p2(libcgen);
    p2.saveSources();
}

void ADDynamics::genLibraryTireFront()
{
    // Front Combined Force
    std::vector<ADCG> x(NX);
    for(int i = 0; i< NX;i++)
        x[i] = 0.0;
    x[3] = 10.;
    Independent(x);
    // dependent variable vector
    std::vector<ADCG> F_comb_f(1);
    // the model equation
    F_comb_f = TireConFront(x);
    ADFun<CGD> fun(x, F_comb_f);

    // generates source code
    ModelCSourceGen<double> cgen(fun, "TireConFront");
    cgen.setCreateJacobian(true);
    ModelLibraryCSourceGen<double> libcgen(cgen);

    // compile source code
    DynamicModelLibraryProcessor<double> p(libcgen, "cppad_cg_TCF");

    GccCompiler<double> compiler;
    std::unique_ptr<DynamicLib<double>> dynamicLibTireConFront = p.createDynamicLibrary(compiler);

    // save to files (not really required)
    SaveFilesModelLibraryProcessor<double> p2(libcgen);
    p2.saveSources();
}

void ADDynamics::genLibraryTireRear()
{
    // Front Combined Force
    std::vector<ADCG> x(NX);
    for(int i = 0; i< NX;i++)
        x[i] = 0.0;
    x[3] = 10.;
    Independent(x);
    // dependent variable vector
    std::vector<ADCG> F_comb_r(1);
    // the model equation
    F_comb_r = TireConRear(x);
    ADFun<CGD> fun(x, F_comb_r);

    // generates source code
    ModelCSourceGen<double> cgen(fun, "TireConRear");
    cgen.setCreateJacobian(true);
    ModelLibraryCSourceGen<double> libcgen(cgen);

    // compile source code
    DynamicModelLibraryProcessor<double> p(libcgen, "cppad_cg_TCR");

    GccCompiler<double> compiler;
    std::unique_ptr<DynamicLib<double>> dynamicLibTireConRear = p.createDynamicLibrary(compiler);

    // save to files (not really required)
    SaveFilesModelLibraryProcessor<double> p2(libcgen);
    p2.saveSources();
}

std::vector<ADCG> ADDynamics::scalerMult(std::vector<ADCG> x,double a)
{
    std::vector<ADCG> y(NX);
    for(int i = 0; i<NX; i++)
        y[i] = a*x[i];
    return y;
}

std::vector<ADCG> ADDynamics::vectorAdd(std::vector<ADCG> x1,std::vector<ADCG> x2)
{
    std::vector<ADCG> y(NX);
    for(int i = 0; i<NX; i++)
        y[i] = x1[i] + x2[i];
    return y;
}

ADCG ADDynamics::getState(std::vector<ADCG> x,int index)
{
    return x[index];
}

ADCG ADDynamics::getInput(std::vector<ADCG> u,int index)
{
    return u[index];
}

ADCG ADDynamics::getSlipAngleFront(std::vector<ADCG> x)
{
    ADCG vx = getState(x,si_index.vx);
    // compute slip angels given current state
    return atan((getState(x,si_index.vy) + getState(x,si_index.r)*param_.lf)/vx)
            - getState(x,si_index.delta);
}

ADCG ADDynamics::getSlipAngleRear(std::vector<ADCG> x)
{
    // compute slip angels given current state
    ADCG vx = getState(x,si_index.vx);
    return atan((getState(x,si_index.vy) - getState(x,si_index.r)*param_.lr)/vx);
}

TireForces ADDynamics::getForceFront(std::vector<ADCG> x)
{
    ADCG alpha_f = getSlipAngleFront(x);
    NormalForces f_normal = getForceNormalDyn(x);

    ADCG F_y = f_normal.F_N_front*param_.Df * sin(param_.Cf * atan(param_.Bf * alpha_f ));
    ADCG F_x = -param_.CBf*getState(x,si_index.B) - param_.Cr0*0.5;

    return {F_x,F_y};
}

TireForces ADDynamics::getForceRear(std::vector<ADCG> x)
{
    ADCG alpha_r = getSlipAngleRear(x);
    NormalForces f_normal  = getForceNormalDyn(x);

    ADCG F_y =  f_normal.F_N_rear*param_.Dr * sin(param_.Cr * atan(param_.Br * alpha_r ));
    ADCG F_x = param_.Cm1*getState(x,si_index.D) - param_.Cm2*getState(x,si_index.D)*getState(x,si_index.vx)
                - param_.CBr*getState(x,si_index.B) - param_.Cr0*0.5;

    return {F_x,F_y};
}

ADCG ADDynamics::getForceFriction(std::vector<ADCG> x)
{
    return -0.5*param_.rho*param_.S*param_.Cr2*getState(x,si_index.vx)*getState(x,si_index.vx);
}

NormalForces ADDynamics::getForceNormalStatic(void)
{
    ADCG F_N_front = (ADCG)(param_.lr/(param_.lf + param_.lr)*param_.m*param_.g);
    ADCG F_N_rear = (ADCG)(param_.lf/(param_.lf + param_.lr)*param_.m*param_.g);
    return {F_N_front,F_N_rear};
}

NormalForces ADDynamics::getForceNormalDyn(std::vector<ADCG> x)
{
    // including aero
    ADCG vx = getState(x,si_index.vx);

    ADCG F_N_front = param_.aero_split_front*
                    0.5*param_.rho*param_.S*param_.Cl*vx*vx
                    + param_.lr/(param_.lf + param_.lr)*param_.m*param_.g;
    ADCG F_N_rear = (1.0 - param_.aero_split_front)*
                    0.5*param_.rho*param_.S*param_.Cl*vx*vx
                    + param_.lf/(param_.lf + param_.lr)*param_.m*param_.g;
    return {F_N_front,F_N_rear};
}

std::vector<ADCG> ADDynamics::dx(std::vector<ADCG> x,std::vector<ADCG> u)
{
    ADCG phi = getState(x,si_index.phi);
    ADCG vx = getState(x,si_index.vx);
    ADCG vy = getState(x,si_index.vy);
    ADCG r  = getState(x,si_index.r);
    ADCG D = getState(x,si_index.D);
    ADCG delta = getState(x,si_index.delta);
    ADCG vs = getState(x,si_index.vs);

    ADCG dD = getInput(u,si_index.dD);
    ADCG dB = getInput(u,si_index.dB);
    ADCG dDelta = getInput(u,si_index.dDelta);
    ADCG dVs = getInput(u,si_index.dVs);

    TireForces tire_forces_front = getForceFront(x);
    TireForces tire_forces_rear = getForceRear(x);
    ADCG friction_force = getForceFriction(x);

    std::vector<ADCG> dx(NX);
    dx[0] = vx*cos(phi) - vy*sin(phi);
    dx[1] = vy*cos(phi) + vx*sin(phi);
    dx[2] = r;
    dx[3] = 1.0/param_.m*(tire_forces_rear.F_x + friction_force
                        - tire_forces_front.F_y*sin(delta) 
                        + tire_forces_front.F_x*cos(delta) 
                        + param_.m*vy*r);
    dx[4] = 1.0/param_.m*(tire_forces_rear.F_y
                        + tire_forces_front.F_x*sin(delta) 
                        + tire_forces_front.F_y*cos(delta) 
                        - param_.m*vx*r);
    dx[5] = 1.0/param_.Iz*(- tire_forces_rear.F_y*param_.lr
                        +(
                            tire_forces_front.F_x*sin(delta) 
                            + tire_forces_front.F_y*cos(delta) 
                        )*param_.lf);
    dx[6] = vs;
    dx[7] = dD;
    dx[8] = dB;
    dx[9] = dDelta;
    dx[10] = dVs;

    return dx;
}

std::vector<ADCG> ADDynamics::RK4(std::vector<ADCG> z)
{

    std::vector<ADCG> state(NX);
    std::vector<ADCG> input(NU);
    std::vector<ADCG> y(NX);

    std::vector<ADCG> k1(NX);
    std::vector<ADCG> k2(NX);
    std::vector<ADCG> k3(NX);
    std::vector<ADCG> k4(NX);

    for(int i = 0;i<NX;i++)
    {
        state[i] = z[i];
    }
    for(int i = 0;i<NU;i++)
    {
        input[i] = z[i+NX];
    }

    k1 = dx(state,input);
    k2 = dx(vectorAdd(state,scalerMult(k1,Ts_*0.5)),input);
    k3 = dx(vectorAdd(state,scalerMult(k2,Ts_*0.5)),input);
    k4 = dx(vectorAdd(state,scalerMult(k3,Ts_)),input);

    for(int i = 0;i<NX;i++)
        y[i] = state[i] + Ts_*(k1[i]/6.0 + k2[i]/3.0 + k3[i]/3.0 + k4[i]/6.0);

    return y;
}
std::vector<ADCG> ADDynamics::TireConFront(std::vector<ADCG> x)
{
    TireForces tire_forces_front = getForceFront(x);
    NormalForces f_normal = getForceNormalStatic();
    NormalForces f_normal_dyn = getForceNormalDyn(x);

    ADCG comb_force = ((param_.e_long*param_.e_long*tire_forces_front.F_x*tire_forces_front.F_x)
                     + (tire_forces_front.F_y*tire_forces_front.F_y)
                     - (param_.e_eps*param_.e_eps*f_normal_dyn.F_N_front*f_normal_dyn.F_N_front*param_.Df*param_.Df))
                     /(f_normal.F_N_front*f_normal.F_N_front);

    return {comb_force};
}

std::vector<ADCG> ADDynamics::TireConRear(std::vector<ADCG> x)
{
    TireForces tire_forces_rear = getForceRear(x);
    NormalForces f_normal = getForceNormalStatic();
    NormalForces f_normal_dyn = getForceNormalDyn(x);

    ADCG comb_force = ((param_.e_long*param_.e_long*tire_forces_rear.F_x*tire_forces_rear.F_x)
                + (tire_forces_rear.F_y*tire_forces_rear.F_y)
                - (param_.e_eps*param_.e_eps*f_normal_dyn.F_N_rear*f_normal_dyn.F_N_rear*param_.Dr*param_.Dr))
                /(f_normal.F_N_rear*f_normal.F_N_rear);

    return {comb_force};
}
std::vector<ADCG> ADDynamics::f_dyn(std::vector<ADCG> z)
{
    std::vector<ADCG> state(NX);
    std::vector<ADCG> input(NU);
    for(int i = 0;i<NX;i++)
    {
        state[i] = z[i];
    }
    for(int i = 0;i<NU;i++)
    {
        input[i] = z[i+NX];
    }

    return dx(state,input);
}

}