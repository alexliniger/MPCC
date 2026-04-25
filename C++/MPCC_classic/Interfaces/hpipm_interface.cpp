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

#include "hpipm_interface.h"
#include "MPC/mpc.h"
namespace mpcc{
void HpipmInterface::setDynamics(std::array<Stage,N+1> &stages,const State &x0)
{
    b0_ = (stages[0].lin_model.A*stateToVector(x0)+stages[0].lin_model.g);
    for(int i=0;i<N;i++)
    {
        if(i==0)
        {
            hA_[i]  = nullptr;
            hB_[i]  = stages[i].lin_model.B.data();
            hb_[i]  = b0_.data();

            nx_[i]  = 0;
            nu_[i]  = NU;
        }
        else
        {
            hA_[i]  = stages[i].lin_model.A.data();
            hB_[i]  = stages[i].lin_model.B.data();
            hb_[i]  = stages[i].lin_model.g.data();

            nx_[i]   = NX;
            nu_[i]   = NU;
        }
    }
    nx_[N]   = NX;
    nu_[N]   = 0;
}

void HpipmInterface::setCost(std::array<Stage,N+1> &stages)
{
    for(int i=0;i<=N;i++) {
        hQ_[i] = stages[i].cost_mat.Q.data();
        hR_[i] = stages[i].cost_mat.R.data();
        hS_[i] = stages[i].cost_mat.S.data();

        hq_[i] = stages[i].cost_mat.q.data();
        hr_[i] = stages[i].cost_mat.r.data();

        if(stages[i].ns != 0)
        {
            hZl_[i] = stages[i].cost_mat.Z.data();
            hZu_[i] = stages[i].cost_mat.Z.data();
            hzl_[i] = stages[i].cost_mat.z.data();
            hzu_[i] = stages[i].cost_mat.z.data();
        }
        else
        {
            hZl_[i] = nullptr;
            hZu_[i] = nullptr;
            hzl_[i] = nullptr;
            hzu_[i] = nullptr;
        }

    }
}

void HpipmInterface::setBounds(std::array<Stage,N+1> &stages,const State &x0)
{

    nbu_[0] = 0;
    hpipm_bounds_[0].idx_u.resize(0);
    hpipm_bounds_[0].lower_bounds_u.resize(0);
    hpipm_bounds_[0].upper_bounds_u.resize(0);
    for(int j=0;j<NU;j++)
    {
        if(stages[0].l_bounds_u(j)>-INF && stages[0].u_bounds_u(j) < INF)
        {
            nbu_[0]++;
            hpipm_bounds_[0].idx_u.push_back(j);
            hpipm_bounds_[0].lower_bounds_u.push_back(stages[0].l_bounds_u(j));
            hpipm_bounds_[0].upper_bounds_u.push_back(stages[0].u_bounds_u(j));
        }
    }
    nbx_[0] = 0;
    hidxbx_[0] = nullptr;
    hidxbu_[0] = hpipm_bounds_[0].idx_u.data();

    hlbx_[0] = nullptr;
    hubx_[0] = nullptr;
    hlbu_[0] = hpipm_bounds_[0].lower_bounds_u.data();
    hubu_[0] = hpipm_bounds_[0].upper_bounds_u.data();

    for(int i=1;i<=N;i++)
    {
        hpipm_bounds_[i].idx_u.resize(0);
        hpipm_bounds_[i].lower_bounds_u.resize(0);
        hpipm_bounds_[i].upper_bounds_u.resize(0);
        nbu_[i] = 0;
        for(int j=0;j<NU;j++)
        {
            if(stages[i].l_bounds_u(j)>-INF && stages[i].u_bounds_u(j) < INF)
            {
                nbu_[i]++;
                hpipm_bounds_[i].idx_u.push_back(j);
                hpipm_bounds_[i].lower_bounds_u.push_back(stages[i].l_bounds_u(j));
                hpipm_bounds_[i].upper_bounds_u.push_back(stages[i].u_bounds_u(j));
            }
        }
        hpipm_bounds_[i].idx_x.resize(0);
        hpipm_bounds_[i].lower_bounds_x.resize(0);
        hpipm_bounds_[i].upper_bounds_x.resize(0);
        nbx_[i] = 0;
        for(int j=0;j<NX;j++)
        {
            if(stages[i].l_bounds_x(j)>-INF && stages[i].u_bounds_x(j) < INF)
            {
                nbx_[i]++;
                hpipm_bounds_[i].idx_x.push_back(j);
                hpipm_bounds_[i].lower_bounds_x.push_back(stages[i].l_bounds_x(j));
                hpipm_bounds_[i].upper_bounds_x.push_back(stages[i].u_bounds_x(j));
            }
        }
        hidxbx_[i] = hpipm_bounds_[i].idx_x.data();
        hidxbu_[i] = hpipm_bounds_[i].idx_u.data();
        hlbx_[i] = hpipm_bounds_[i].lower_bounds_x.data();
        hubx_[i] = hpipm_bounds_[i].upper_bounds_x.data();
        hlbu_[i] = hpipm_bounds_[i].lower_bounds_u.data();
        hubu_[i] = hpipm_bounds_[i].upper_bounds_u.data();

    }

    nbu_[N] = 0;
    hidxbu_[N] = nullptr;
    hlbu_[N] = nullptr;
    hubu_[N] = nullptr;
}

void HpipmInterface::setPolytopicConstraints(std::array<Stage,N+1> &stages)
{
    for(int i=0;i<=N;i++) {
        ng_[i] = stages[i].ng;
        if (stages[i].ng > 0) {
            hC_[i] = stages[i].constrains_mat.C.data();
            hD_[i] = stages[i].constrains_mat.D.data();

            hlg_[i] = stages[i].constrains_mat.dl.data();
            hug_[i] = stages[i].constrains_mat.du.data();
        }
        else
        {
            hC_[i] = nullptr;
            hD_[i] = nullptr;

            hlg_[i] = nullptr;
            hug_[i] = nullptr;
        }
    }
}

void HpipmInterface::setSoftConstraints(std::array<Stage,N+1> &stages)
{
    for(int i=0;i<=N;i++)
    {
        hpipm_bounds_[i].idx_s.resize(0);
        if (stages[i].ns != 0)
        {
            nsbx_[i] = 0;
            nsbu_[i] = 0;
            nsg_[i] = stages[i].ns;

            for(int j=0;j<stages[i].ns;j++)
            {
                hpipm_bounds_[i].idx_s.push_back(j+nbx_[i]+nbu_[i]);
            }

            hidxs_[i] = hpipm_bounds_[i].idx_s.data();
            hlls_[i] = stages[i].l_bounds_s.data();
            hlus_[i] = stages[i].u_bounds_s.data();
        }
        else
        {
            nsbx_[i] = 0;
            nsbu_[i] = 0;
            nsg_[i] = 0;
            hidxs_[i] = nullptr;
            hlls_[i] = nullptr;
            hlus_[i] = nullptr;
        }
    }
}


std::array<OptVariables,N+1> HpipmInterface::solveMPC(std::array<Stage,N+1> &stages, const State &x0, int *status)
{
    setDynamics(stages,x0);
    setCost(stages);
    setBounds(stages,x0);
    setPolytopicConstraints(stages);
    setSoftConstraints(stages);
//    print_data();

    std::array<OptVariables,N+1> opt_solution = Solve(status);
    opt_solution[0].xk = x0;

    return opt_solution;
}

std::array<OptVariables,N+1> HpipmInterface::Solve(int *status)
{
    // ocp qp dim
    int dim_size = d_ocp_qp_dim_memsize(N);
    void *dim_mem = malloc(dim_size);

    struct d_ocp_qp_dim dim;
    d_ocp_qp_dim_create(N, &dim, dim_mem);

    d_ocp_qp_dim_set_all(nx_, nu_, nbx_, nbu_, ng_, nsbx_, nsbu_, nsg_, &dim);
    // ocp qp
    int qp_size = d_ocp_qp_memsize(&dim);
    void *qp_mem = malloc(qp_size);

    struct d_ocp_qp qp;
    d_ocp_qp_create(&dim, &qp, qp_mem);
    d_ocp_qp_set_all(hA_, hB_, hb_, hQ_, hS_, hR_, hq_, hr_,
                     hidxbx_, hlbx_, hubx_, hidxbu_, hlbu_, hubu_,
                     hC_, hD_, hlg_, hug_, hZl_, hZu_, hzl_, hzu_,
                     hidxs_, hlls_, hlus_, &qp);

    // ocp qp sol
    int qp_sol_size = d_ocp_qp_sol_memsize(&dim);
    void *qp_sol_mem = malloc(qp_sol_size);

    struct d_ocp_qp_sol qp_sol;
    d_ocp_qp_sol_create(&dim, &qp_sol, qp_sol_mem);

    // ipm arg

    int ipm_arg_size = d_ocp_qp_ipm_arg_memsize(&dim);
    printf("\nipm arg size = %d\n", ipm_arg_size);
    void *ipm_arg_mem = malloc(ipm_arg_size);

    struct d_ocp_qp_ipm_arg arg;
    d_ocp_qp_ipm_arg_create(&dim, &arg, ipm_arg_mem);

//    enum hpipm_mode mode = SPEED_ABS;
    enum hpipm_mode mode = SPEED;
//    enum hpipm_mode mode = BALANCE;
//    enum hpipm_mode mode = ROBUST;

//    int mode = 1;
    double mu0 = 1e2;
    int iter_max = 30;
    double tol_stat = 1e-6;
    double tol_eq = 1e-6;
    double tol_ineq = 1e-6;
    double tol_comp = 1e-6;
    double reg_prim = 1e-12;
    int warm_start = 0;
    int pred_corr = 1;
    int ric_alg = 0;

    d_ocp_qp_ipm_arg_set_default(mode, &arg);

    // d_ocp_qp_ipm_arg_set_mu0(&mu0, &arg);
    d_ocp_qp_ipm_arg_set_iter_max(&iter_max, &arg);
//    d_ocp_qp_ipm_arg_set_tol_stat(&tol_stat, &arg);
//    d_ocp_qp_ipm_arg_set_tol_eq(&tol_eq, &arg);
//    d_ocp_qp_ipm_arg_set_tol_ineq(&tol_ineq, &arg);
//    d_ocp_qp_ipm_arg_set_tol_comp(&tol_comp, &arg);
//    d_ocp_qp_ipm_arg_set_reg_prim(&reg_prim, &arg);
//    d_ocp_qp_ipm_arg_set_warm_start(&warm_start, &arg);
//    d_ocp_qp_ipm_arg_set_pred_corr(&pred_corr, &arg);
//    d_ocp_qp_ipm_arg_set_ric_alg(&ric_alg, &arg);

    // ipm

    int ipm_size = d_ocp_qp_ipm_ws_memsize(&dim, &arg);
//    printf("\nipm size = %d\n", ipm_size);
    void *ipm_mem = malloc(ipm_size);

    struct d_ocp_qp_ipm_ws workspace;
    d_ocp_qp_ipm_ws_create(&dim, &arg, &workspace, ipm_mem);

    int hpipm_return; // 0 normal; 1 max iter; 2 linesearch issues?

    struct timeval tv0, tv1;

    gettimeofday(&tv0, nullptr); // start
    d_ocp_qp_ipm_solve(&qp, &qp_sol, &arg, &workspace);
    d_ocp_qp_ipm_get_status(&workspace, &hpipm_return);
    gettimeofday(&tv1, nullptr); // stop
    double time_ocp_ipm = (tv1.tv_usec-tv0.tv_usec)/(1e6);

    printf("comp time = %f\n", time_ocp_ipm);
    printf("exitflag %d\n", hpipm_return);
    printf("ipm iter = %d\n", workspace.iter);

    // extract and print solution
    int ii;
    int nu_max = nu_[0];
    for(ii=1; ii<=N; ii++)
        if(nu_[ii]>nu_max)
            nu_max = nu_[ii];
    double *u = (double*)malloc(nu_max*sizeof(double));
//    printf("\nu = \n");
    for(ii=0; ii<=N; ii++)
    {
        d_ocp_qp_sol_get_u(ii, &qp_sol, u);
//        d_print_mat(1, nu_[ii], u, 1);
    }
    int nx_max = nx_[0];
    for(ii=1; ii<=N; ii++)
        if(nx_[ii]>nx_max)
            nx_max = nx_[ii];
    double *x = (double*)malloc(nx_max*sizeof(double));
//    printf("\nx = \n");
    for(ii=0; ii<=N; ii++)
    {
        d_ocp_qp_sol_get_x(ii, &qp_sol, x);
//        d_print_mat(1, nx_[ii], x, 1);
    }


    std::array<OptVariables,N+1> optimal_solution;
    optimal_solution[0].xk.setZero();
    for(int i=1;i<=N;i++){
        d_ocp_qp_sol_get_x(i, &qp_sol, x);
        optimal_solution[i].xk = arrayToState(x);
    }

    for(int i=0;i<N;i++){
        d_ocp_qp_sol_get_u(i, &qp_sol, u);
        optimal_solution[i].uk = arrayToInput(u);
    }
    optimal_solution[N].uk.setZero();

    /************************************************
    * print ipm statistics
    ************************************************/

    int iter; d_ocp_qp_ipm_get_iter(&workspace, &iter);
    double res_stat; d_ocp_qp_ipm_get_max_res_stat(&workspace, &res_stat);
    double res_eq; d_ocp_qp_ipm_get_max_res_eq(&workspace, &res_eq);
    double res_ineq; d_ocp_qp_ipm_get_max_res_ineq(&workspace, &res_ineq);
    double res_comp; d_ocp_qp_ipm_get_max_res_comp(&workspace, &res_comp);
    double *stat; d_ocp_qp_ipm_get_stat(&workspace, &stat);
    int stat_m; d_ocp_qp_ipm_get_stat_m(&workspace, &stat_m);

//    printf("\nipm return = %d\n", hpipm_return);
//    printf("\nipm residuals max: res_g = %e, res_b = %e, res_d = %e, res_m = %e\n", res_stat, res_eq, res_ineq, res_comp);
//
//    printf("\nipm iter = %d\n", iter);
//    printf("\nalpha_aff\tmu_aff\t\tsigma\t\talpha\t\tmu\t\tres_stat\tres_eq\t\tres_ineq\tres_comp\n");
//    d_print_exp_tran_mat(stat_m, iter+1, stat, stat_m);

    free(dim_mem);
    free(qp_mem);
    free(qp_sol_mem);
    free(ipm_arg_mem);
    free(ipm_mem);

    free(u);
    free(x);

    *status = hpipm_return;

    return optimal_solution;
}

void HpipmInterface::print_data(){

    for(int k=0; k<=N ; k++) {
        if (k != N) {
            std::cout << "A_"<< k << " = " << std::endl;
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NX; j++) {
                    std::cout << hA_[k][i + j * NX] << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "B_"<< k << " = " << std::endl;
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < NU; j++) {
                    std::cout << hB_[k][i + j * NX] << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "b_" << k << " = " << std::endl;
            for (int i = 0; i < NX; i++) {
                for (int j = 0; j < 1; j++) {
                    std::cout << hb_[k][i + j * NX] << " ";
                }
                std::cout << std::endl;
            }
        }

        std::cout << "Q_" << k << " = " << std::endl;
        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NX; j++) {
                std::cout << hQ_[k][i + j * NX] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "R_" << k << " = " << std::endl;
        for (int i = 0; i < NU; i++) {
            for (int j = 0; j < NU; j++) {
                std::cout << hR_[k][i + j * NU] << " ";
            }
            std::cout << std::endl;
        }

        std::cout << "q_" << k << " = " << std::endl;
        for (int i = 0; i < NX; i++) {
            std::cout << hq_[k][i] << " ";
            std::cout << std::endl;
        }

        std::cout << "r_" << k << " = " << std::endl;
        for (int i = 0; i < NU; i++) {
            std::cout << hr_[k][i] << " ";
            std::cout << std::endl;
        }

        std::cout << "sizes" << std::endl;
        std::cout << nx_[k] << std::endl;
        std::cout << nu_[k] << std::endl;
        std::cout << nbx_[k] << std::endl;
        std::cout << nbu_[k] << std::endl;
        std::cout << ng_[k] << std::endl;
    }

}
}