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

#include "mpc.h"

namespace mpcc{
void MPC::setMPCProblem()
{
    for(int i=0;i<=N;i++)
    {
        setStage(initial_guess_[i].xk,initial_guess_[i].uk,i);
    }
}

void MPC::setStage(const State &xk, const Input &uk, const int time_step)
{
    stages_[time_step].nx = NX;
    stages_[time_step].nu = NU;

    if(time_step == 0)
    {
        stages_[time_step].ng = 0;
        stages_[time_step].ns = 0;
    }
    else
    {
        stages_[time_step].ng = 3;
        stages_[time_step].ns = 3;
    }

    State xk_nz = xk;
    xk_nz.vxNonZero();

    stages_[time_step].cost_mat = cost_.getCost(track_,xk_nz,time_step);
    stages_[time_step].lin_model = model_.getLinModel(xk_nz,uk);
    stages_[time_step].constrains_mat = constraints_.getConstraints(track_,xk_nz,uk);
    //TODO better bounds interface
    stages_[time_step].l_bounds_x = bounds_.getBoundsLX();
    stages_[time_step].u_bounds_x = bounds_.getBoundsUX();
    stages_[time_step].l_bounds_u = bounds_.getBoundsLU();
    stages_[time_step].u_bounds_u = bounds_.getBoundsUU();
    stages_[time_step].l_bounds_s = bounds_.getBoundsLS();
    stages_[time_step].u_bounds_s = bounds_.getBoundsUS();

    stages_[time_step].l_bounds_x(si_index.s) = initial_guess_[time_step].xk.s - param.s_trust_region;//*initial_guess_[time_step].xk.vs;
    stages_[time_step].u_bounds_x(si_index.s) = initial_guess_[time_step].xk.s + param.s_trust_region;//*initial_guess_[time_step].xk.vs;
}

void MPC::updateInitialGuess(const State &x0)
{
    for(int i=1;i<N;i++)
        initial_guess_[i-1] = initial_guess_[i];

    initial_guess_[0].xk = x0;
    initial_guess_[0].uk.setZero();

    initial_guess_[N-1].xk = initial_guess_[N-2].xk;
    initial_guess_[N-1].uk.setZero();// = initial_guess_[N-2].uk;

    initial_guess_[N].xk = integrator_.RK4(initial_guess_[N-1].xk,initial_guess_[N-1].uk,TS);
    initial_guess_[N].uk.setZero();

    unwrapInitialGuess();
}

// alternatively OptVariables MPC::unwrapInitialGuess(const OptVariables &initial_guess)
void MPC::unwrapInitialGuess()
{
    double L = track_.getLength();
    for(int i=1;i<=N;i++)
    {
        if((initial_guess_[i].xk.phi - initial_guess_[i-1].xk.phi) < M_PI)
        {
            initial_guess_[i].xk.phi += 2.*M_PI;
        }
        if((initial_guess_[i].xk.phi - initial_guess_[i-1].xk.phi) > M_PI)
        {
            initial_guess_[i].xk.phi -= 2.*M_PI;
        }

        if((initial_guess_[i].xk.s - initial_guess_[i-1].xk.s) > L/2.)
        {
            initial_guess_[i].xk.s -= L;
        }
    }

}

void MPC::generateNewInitialGuess(const State &x0)
{
    initial_guess_[0].xk = x0;
    initial_guess_[0].uk.setZero();

    for(int i = 1;i<=N;i++)
    {
        initial_guess_[i].xk.setZero();
        initial_guess_[i].uk.setZero();

        initial_guess_[i].xk.s = initial_guess_[i-1].xk.s + TS*param.initial_velocity;
        Eigen::Vector2d track_pos_i = track_.getPostion(initial_guess_[i].xk.s);
        Eigen::Vector2d track_dpos_i = track_.getDerivative(initial_guess_[i].xk.s);
        initial_guess_[i].xk.X = track_pos_i(0);
        initial_guess_[i].xk.Y = track_pos_i(1);
        initial_guess_[i].xk.phi = atan2(track_dpos_i(1),track_dpos_i(0));
        initial_guess_[i].xk.vx = param.initial_velocity;
        initial_guess_[i].xk.vs = param.initial_velocity;
    }
    unwrapInitialGuess();
    valid_initial_guess_ = true;
}

std::array<OptVariables,N+1> MPC::sqpSolutionUpdate(const std::array<OptVariables,N+1> &last_solution,
                                                    const std::array<OptVariables,N+1> &current_solution)
{
    //TODO use line search and merit function
    std::array<OptVariables,N+1> updated_solution;
    StateVector updated_x_vec;
    InputVector updated_u_vec;
    for(int i = 0;i<=N;i++)
    {
        updated_x_vec = param.sqp_mixing*stateToVector(current_solution[i].xk)
                        +(1.0-param.sqp_mixing)*stateToVector(last_solution[i].xk);
        updated_u_vec = param.sqp_mixing*inputToVector(current_solution[i].uk)
                        +(1.0-param.sqp_mixing)*inputToVector(last_solution[i].uk);

        updated_solution[i].xk = vectorToState(updated_x_vec);
        updated_solution[i].uk = vectorToInput(updated_u_vec);
    }

    return updated_solution;
}

MPCReturn MPC::runMPC(State &x0)
{
    struct timeval tv0, tv1;
    gettimeofday(&tv0, nullptr); // start
    // TODO: implemment projection onto centerline for x0
    x0.unwrap(track_.getLength());
    if(valid_initial_guess_)
        updateInitialGuess(x0);
    else
        generateNewInitialGuess(x0);

    //TODO: Detect if not solved to optimility for to long and restart with fixed velocity initial guess
    // change valid_initial_guess_ = 0 to restart
    for(int i=0;i<n_sqp_;i++)
    {
        setMPCProblem();
        optimal_solution_ = solver_interface_->solveMPC(stages_,x0);
        initial_guess_ = sqpSolutionUpdate(initial_guess_,optimal_solution_);

    }
    gettimeofday(&tv1, nullptr); // stop
    double time_nmpc = (tv1.tv_usec-tv0.tv_usec)/(1e6);

    return {initial_guess_[0].uk,initial_guess_,time_nmpc};
}

void MPC::setTrack(const Eigen::VectorXd &X, const Eigen::VectorXd &Y){
    track_.gen2DSpline(X,Y);
}

MPC::MPC(int n_sqp)
:valid_initial_guess_(false),solver_interface_(new HpipmInterface())
{
    n_sqp_ = n_sqp;
}
}