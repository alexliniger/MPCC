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

function [ X,U,dU,info ] = YalmipInterface(stage,MPC_vars,ModelParams)
nx = ModelParams.nx;
nu = ModelParams.nu;
tic
x    = sdpvar(nx+nu, MPC_vars.N+1);          % states + initial state; fifth initial state for discretization
u    = sdpvar(nu, MPC_vars.N);           % (front) steering angle
objective = 0;
constraints = [x(:,1) == blkdiag(MPC_vars.Tx,MPC_vars.Tu)*[stage(1).x0;stage(1).u0]];   % initialize initial state

for i = 1:MPC_vars.N
    constraints = [constraints;
                   x(:, i+1) == stage(i).Ak*x(:,i) + stage(i).Bk*u(:,i) + stage(i).gk ; %dynamics
                   stage(i).lb <= [x(:,i);u(:,i)] <= stage(i).ub; % bounds
                   stage(i).lg <= stage(i).Ck*x(:,i) <= stage(i).ug]; %track constraints

    objective = objective + 0.5*(x(:,i)'*stage(i).Qk*x(:,i) + u(:,i)'*stage(i).Rk*u(:,i)) ...
                          + stage(i).fk'*x(:,i);  % cost
 
end
i = MPC_vars.N+1;
objective = objective + 0.5*(x(:,i)'*stage(i).Qk*x(:,i)) + stage(i).fk'*x(:,i);
constraints = [constraints;
               stage(i).lb(1:nx+nu) <= x(:,i) <= stage(i).ub(1:nx+nu); %bounds
               stage(i).lg <= stage(i).Ck*x(:,i) <= stage(i).ug];  % track constraints
yalmipTime = toc           
ops = sdpsettings('solver','ecos','verbose',2);

% solve QP
tic;
exitflag = solvesdp(constraints,objective,ops);
QPtime = toc;

x_opt = double(x);
u_opt = double(u);

% rescale outputs
X= MPC_vars.invTx*x_opt(1:nx,:);
U = MPC_vars.invTu*x_opt(nx+1:end,2:end);
dU = u_opt;

info.exitflag = exitflag;
info.QPtime = QPtime;

end

