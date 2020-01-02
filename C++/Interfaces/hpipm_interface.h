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

#ifndef MPCC_HPIPM_INTERFACE_H
#define MPCC_HPIPM_INTERFACE_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#include <blasfeo_d_aux_ext_dep.h>

#include "hpipm_d_ocp_qp_ipm.h"
#include "hpipm_d_ocp_qp_dim.h"
#include "hpipm_d_ocp_qp.h"
#include "hpipm_d_ocp_qp_sol.h"
#include "hpipm_timing.h"

#include "config.h"
#include "types.h"
#include "Model/model.h"
#include "Cost/cost.h"
#include "Constraints/constraints.h"
#include "Constraints/bounds.h"
#include "solver_interface.h"

#include <array>
#include <vector>

namespace mpcc{
struct OptVariables;
struct Stage;

struct HpipmBound {
    std::vector<int> idx_u;
    std::vector<int> idx_x;
    std::vector<int> idx_s;
    std::vector<double> lower_bounds_u;
    std::vector<double> upper_bounds_u;
    std::vector<double> lower_bounds_x;
    std::vector<double> upper_bounds_x;
};

class HpipmInterface : public SolverInterface {
public:
    std::array<OptVariables,N+1> solveMPC(std::array<Stage,N+1> &stages,const State &x0, int *status);

    ~HpipmInterface(){
        std::cout << "Deleting Hpipm Interface" << std::endl;
    }
private:
    int nx_[N+1];//    -> number of states
    int nu_[N+1];//    -> number of inputs
    int nbx_[N+1];//   -> number of bounds on x
    int nbu_[N+1];//   -> number of bounds on u
    int ng_[N+1];//    -> number of polytopic constratins
    int nsbx_[N+1];//   -> number of slacks variables on x
    int nsbu_[N+1];//   -> number of slacks variables on u
    int nsg_[N+1];//   -> number of slacks variables on polytopic constraints

    // LTV dynamics
    // x_k+1 = A_k x_k + B_k u_k + b_k
    double *hA_[N]; //hA[k] = A_k
    double *hB_[N]; //hB[k] = B_k
    double *hb_[N]; //hb[k] = b_k

    // Cost (without soft constraints)
    // min_x,u sum_k=0^N 1/2*[x_k;u_k]^T*[Q_k , S_k; S_k^T , R_k]*[x_k;u_k] + [q_k; r_k]^T*[x_k;u_k]
    double *hQ_[N+1]; //hQ[k] = Q_k
    double *hS_[N+1]; //hS[k] = S_k
    double *hR_[N+1]; //hR[k] = R_k
    double *hq_[N+1]; //hq[k] = q_k
    double *hr_[N+1]; //hr[k] = r_k

    // Polytopic constraints
    // g_lower,k <= D_k*x_k + C_k*u_k
    // D_k*x_k + C_k*u_k  <= g_upper,k
    double *hlg_[N+1]; //hlg[k] =  g_lower,k
    double *hug_[N+1]; //hug[k] =  g_upper,k
    double *hC_[N+1]; //hC[k] = C_k
    double *hD_[N+1]; //hD[k] = D_k

    // General bounds
    // x_lower,k <= x_k <= x_upper,k
    // hidxbx can be used to select bounds on a subset of states
    int *hidxbx_[N+1]; // hidxbx[k] = {0,1,2,...,nx} for bounds on all inputs and states
    double *hlbx_[N+1]; // x_lower,k
    double *hubx_[N+1]; //x_upper,k
    // u_lower,k <= u_k <=  u_upper,k
    // hidxbu can be used to select bounds on a subset of inputs
    int *hidxbu_[N+1]; // hidxbuk] = {0,1,2,...,nu} for bounds on all inputs and states
    double *hlbu_[N+1]; // u_lower,k
    double *hubu_[N+1]; // u_upper,k

    // Cost (only soft constriants)
    // s_lower,k -> slack variable of lower polytopic constraint (3) + lower bounds
    // s_upper,k -> slack variable of upper polytopic constraint (4) + upper bounds
    // min_x,u sum_k=0^N 1/2*[s_lower,k;s_upper,k]^T*[Z_lower,k , 0; 0 , Z_upper,k]*[s_lower,k;s_upper,k] + [z_lower,k; z_upper,k]^T*[s_lower,k;s_upper,k]
    double *hZl_[N+1]; // hZl[k] = Z_lower,k
    double *hZu_[N+1]; // hZu[k] = Z_upper,k
    double *hzl_[N+1]; // hzl[k] = z_lower,k
    double *hzu_[N+1]; // hzu[k] = z_upper,k

    // Bounds on the soft constraint multipliers
    // s_lower,k >= s_lower,bound,k
    // s_upper,k >= s_upper,bound,k
    double *hlls_[N+1];
    double *hlus_[N+1];
    // index of the bounds and constraints that are softened
    // order is not really clear
    int *hidxs_[N+1];

    //bounds that are different to stages bounds and need to be stored somewhere such the a pointer can point
    std::array<HpipmBound,N+1> hpipm_bounds_;
    Eigen::Matrix<double,NX,1> b0_;

    void setDynamics(std::array<Stage,N+1> &stages,const State &x0);
    void setCost(std::array<Stage,N+1> &stages);
    void setBounds(std::array<Stage,N+1> &stages,const State &x0);
    void setPolytopicConstraints(std::array<Stage,N+1> &stages);
    void setSoftConstraints(std::array<Stage,N+1> &stages);

    std::array<OptVariables,N+1> Solve(int *status);

    void print_data();
};
}
#endif //MPCC_HPIPM_INTERFACE_H