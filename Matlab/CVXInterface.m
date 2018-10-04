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

function [ X,U,dU,info ] = CVXInterface(stage,MPC_vars,ModelParams)
nx = ModelParams.nx;
nu = ModelParams.nu;
tic

cvx_begin
    variable x((nx+nu)*(MPC_vars.N+1));
    variable u(nu*MPC_vars.N);
    
    objective = 0;

    for i = 1:MPC_vars.N
        objective = objective + 0.5*(quad_form(x((i-1)*(nx+nu)+[1:nx+nu]),stage(i).Qk) + quad_form(u((i-1)*nu+[1:nu]),stage(i).Rk)) + stage(i).fk'*x((i-1)*(nx+nu)+[1:nx+nu]);  % cost
    end
    i = MPC_vars.N+1;
    objective = objective + 0.5*(quad_form(x((i-1)*(nx+nu)+[1:nx+nu]),stage(i).Qk)) + stage(i).fk'*x((i-1)*(nx+nu)+[1:nx+nu]);

    minimize( objective )

    subject to              
        x(1:nx+nu) == blkdiag(MPC_vars.Tx,MPC_vars.Tu)*[stage(1).x0;stage(1).u0];   % initialize initial state

        for i = 1:MPC_vars.N
           x((i)*(nx+nu)+[1:nx+nu]) == stage(i).Ak*x((i-1)*(nx+nu)+[1:nx+nu]) + stage(i).Bk*u((i-1)*nu+[1:nu]) + stage(i).gk ; %dynamics
           stage(i).lb <= [x((i-1)*(nx+nu)+[1:nx+nu]);u((i-1)*nu+[1:nu])] <= stage(i).ub; % bounds
           stage(i).lg <= stage(i).Ck*x((i-1)*(nx+nu)+[1:nx+nu]) <= stage(i).ug; %track constraints

        end
        i = MPC_vars.N+1;
        stage(i).lb(1:nx+nu) <= x((i-1)*(nx+nu)+[1:nx+nu]) <= stage(i).ub(1:nx+nu); %bounds
        stage(i).lg <= stage(i).Ck*x((i-1)*(nx+nu)+[1:nx+nu]) <= stage(i).ug;  % track constraints
    
cvx_end   

QPtime = toc();

x_opt = reshape(x,nx+nu,MPC_vars.N+1);
u_opt = reshape(u,nu,MPC_vars.N);

% rescale outputs
X= MPC_vars.invTx*x_opt(1:nx,:);
U = MPC_vars.invTu*x_opt(nx+1:end,2:end);
dU = u_opt;

if cvx_status == "Solved"
    info.exitflag = 0;
else
    info.exitflag = 1;
end
info.QPtime = QPtime;

cvx_clear

end

