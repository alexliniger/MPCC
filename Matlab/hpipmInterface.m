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


% check that env.sh has been run
env_run = getenv('ENV_RUN');
if (~strcmp(env_run, 'true'))
	disp('ERROR: env.sh has not been sourced! Before using the HPIPM solver, run in the shell:');
	disp('source env.sh');
	disp('and then launch matlab or octave from the same shell.');
	disp('(tested on linux, untested on windows as of now).');
	return;
end


nx = ModelParams.nx;
nu = ModelParams.nu;
qpTotal = tic;

%import hpipm_matlab.*

% dims
N = MPC_vars.N;

dims = hpipm_ocp_qp_dim(N);

dims.set('nx', (nx+nu), 0, N);
dims.set('nu', nu, 0, N-1);

dims.set('nbx', (nx+nu), 0, N);
dims.set('nbu', nu, 0, N-1);

dims.set('ng', 0, 0);
dims.set('ng', 1, 1, N);

%dims.print_C_struct();


% qp
qp = hpipm_ocp_qp(dims);
%% Equility Constraints
x0 = blkdiag(MPC_vars.Tx,MPC_vars.Tu)*[stage(1).x0;stage(1).u0];
for i = 0:N-1
   qp.set('A', stage(i+1).Ak, i); 
   qp.set('B', stage(i+1).Bk, i); 
   qp.set('b', stage(i+1).gk, i); 
end

%% Cost
for i = 0:N
    qp.set('Q', stage(i+1).Qk, i);
    qp.set('q', stage(i+1).fk, i);
    if i<N
        qp.set('R', stage(i+1).Rk, i); 
    end
end
%% Constraints
for i = 1:N
    qp.set('C', stage(i+1).Ck, i);
    qp.set('lg', stage(i+1).lg, i); 
    qp.set('ug', stage(i+1).ug, i); 
end
%% Bounds
%qp.print_C_struct();
for i = 0:N
    qp.set('Jbx', eye(nx+nu), i)
    if i == 0
        qp.set('lbx', x0, 0)
        qp.set('ubx', x0, 0)
    else
        qp.set('lbx', stage(i+1).lb(1:nx+nu), i)
        qp.set('ubx', stage(i+1).ub(1:nx+nu), i)
    end
    
    if i<N
        qp.set('Jbu', eye(nu), i)
        qp.set('lbu', stage(i+1).lb(nx+nu+1:nx+nu+nu), i)
        qp.set('ubu', stage(i+1).ub(nx+nu+1:nx+nu+nu), i)
    end
end
    
%qp.print_C_struct();


%%
qp_sol = hpipm_ocp_qp_sol(dims);


%% set up solver arg
%mode = 'speed_abs';
mode = 'speed';
%mode = 'balance';
%mode = 'robust';
% create and set default arg based on mode
arg = hpipm_ocp_qp_solver_arg(dims, mode);

arg.set('mu0', 1e0);
arg.set('iter_max', 200);
arg.set('tol_stat', 1e-6);
arg.set('tol_eq', 1e-6);
arg.set('tol_ineq', 1e-6);
arg.set('tol_comp', 1e-5);
arg.set('reg_prim', 1e-12);


% set up solver
solver = hpipm_ocp_qp_solver(dims, arg);


% solve qp
qptime = tic;
solver.solve(qp, qp_sol);
tmp_time = toc(qptime);

return_flag = solver.get('status');

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
    u_opt(:,i+1) = qp_sol.get('u', i);
end
for i=0:N
	x_opt(:,i+1) = qp_sol.get('x', i);
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



if is_octave()
	% directly call destructor for octave 4.2.2 (ubuntu 18.04) + others ???
	if strcmp(version(), '4.2.2')
		delete(dims);
		delete(qp);
		delete(qp_sol);
		delete(arg);
		delete(solver);
	end
end



end

