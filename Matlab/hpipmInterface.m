% Copyright (C) 2018, ETH Zurich, D-ITET, Kenneth Kuchera, Alexander Liniger
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ X,U,dU,info ] = hpipmInterface(stage,MPC_vars,ModelParams)
nx = ModelParams.nx;
nu = ModelParams.nu;
qpTotal = tic;

import hpipm_matlab.*

% dims
N = MPC_vars.N;

dims = hpipm_ocp_qp_dim(N);

dims.set_nx((nx+nu)*ones(1,N+1));
dims.set_nu(nu*ones(1,N));

dims.set_nbx((nx+nu)*ones(1,N+1));
dims.set_nbu(nu*ones(1,N));

dims.set_ng([0,1*ones(1,N)]);

% qp
qp = hpipm_ocp_qp(dims);
%% Equility Constraints
x0 = blkdiag(MPC_vars.Tx,MPC_vars.Tu)*[stage(1).x0;stage(1).u0];
for i = 0:N-1
   qp.set_A(stage(i+1).Ak,i); 
   qp.set_B(stage(i+1).Bk,i); 
   qp.set_b(stage(i+1).gk,i); 
end

%% Cost
for i = 0:N
    qp.set_Q(stage(i+1).Qk,i);
    qp.set_q(stage(i+1).fk,i);
    if i<N
        qp.set_R(stage(i+1).Rk,i); 
    end
end
%% Constraints
for i = 1:N
    qp.set_C(stage(i+1).Ck,i);
    qp.set_lg(stage(i+1).lg,i); 
    qp.set_ug(stage(i+1).ug,i); 
end
%% Bounds
for i = 0:N
    qp.set_Jx(eye(nx+nu),i)
    if i == 0
        qp.set_lx(x0,0)
        qp.set_ux(x0,0)
    else
        qp.set_lx(stage(i+1).lb(1:nx+nu),i)
        qp.set_ux(stage(i+1).ub(1:nx+nu),i)
    end
    
    if i<N
        qp.set_Ju(eye(nu),i)
        qp.set_lu(stage(i+1).lb(nx+nu+1:end),i)
        qp.set_uu(stage(i+1).ub(nx+nu+1:end),i)
    end
end
    
%%
qp_sol = hpipm_ocp_qp_sol(dims);

% set up solver arg
arg = hpipm_ocp_qp_solver_arg(dims);

arg.set_mu0(1e0);
arg.set_iter_max(200);
arg.set_tol_stat(1e-6);
arg.set_tol_eq(1e-6);
arg.set_tol_ineq(1e-6);
arg.set_tol_comp(1e-5);
arg.set_reg_prim(1e-12);

% set up solver
solver = hpipm_ocp_qp_solver(dims, arg);

% solve qp
qptime = tic;
return_flag = solver.solve(qp, qp_sol);
tmp_time = toc(qptime);
fprintf('solve time %e\n', tmp_time);

fprintf('HPIPM returned with flag %d ', return_flag);

if return_flag==0
    fprintf('-> QP solved\n')
%     qp_sol.print_C_struct()
else
    fprintf('-> Solver failed!\n')
end


% extract and print sol
u_opt = zeros(nu,N);
x_opt = zeros(nx+nu,N);
for i=0:N-1
    u_opt(:,i+1) = qp_sol.get_u(i);
end
for i=0:N
	x_opt(:,i+1) = qp_sol.get_x(i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
QPtime = toc(qpTotal);

% rescale outputs
X = MPC_vars.invTx*x_opt(1:nx,:);
U = MPC_vars.invTu*x_opt(nx+1:end,2:end);
dU = u_opt;

if return_flag == 0
    info.exitflag = 0;
else
    info.exitflag = 1;
end
info.QPtime = tmp_time;


end

