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

function [xTemp,uTemp] = augState(x,u,x0,MPC_vars,ModelParams,tl)
    
    nx = ModelParams.nx;
    nu = ModelParams.nu;
    N = MPC_vars.N;
    Ts = MPC_vars.Ts;
    indexPhi = ModelParams.stateindex_phi;
    indexTheta = ModelParams.stateindex_theta;

    xTemp = zeros(nx,N+1);
    uTemp = zeros(nu,N);
    
    xTemp(:,1) = x0;
    uTemp(:,1) = u(:,2);
    for j=2:N-1
        xTemp(:,j) = x(:,j+1);
        uTemp(:,j) = u(:,j+1);
    end
    j = N;
    xTemp(:,j) = x(:,j+1);
    uTemp(:,j) = u(:,j);
    
    j = N+1;
    xTemp(:,j) = SimTimeStep(x(:,N+1),u(:,N),Ts,ModelParams);
    
    if xTemp(indexPhi,1) - xTemp(indexPhi,2) > pi
        xTemp(indexPhi,2:end) = xTemp(indexPhi,2:end)+2*pi;
    end
    if xTemp(indexPhi,1) - xTemp(indexPhi,2) < -pi
        xTemp(indexPhi,2:end) = xTemp(indexPhi,2:end)-2*pi;
    end
    
    if xTemp(indexTheta,1) - xTemp(indexTheta,2) < -0.75*tl
        xTemp(indexTheta,2:end) = xTemp(indexTheta,2:end)-tl;
    end
    
end