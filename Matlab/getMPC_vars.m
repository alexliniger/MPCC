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

function MPC_vars = getMPC_vars(CarModel)


if CarModel == "ORCA"
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MPC settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prediction horizon
    MPC_vars.N = 40;
    % sampling time
    MPC_vars.Ts = 0.02;
    % used model (TODO incorparate new models)
    MPC_vars.ModelNo = 1;
    % use bounds on all opt variables (TODO implement selective bounds)
    MPC_vars.fullBound = 1;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % state-input scaling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % normalization matricies (scale states and inputs to ~ +/- 1
    MPC_vars.Tx = diag(1./[3,3,2*pi,4,2,7,30]);
    MPC_vars.Tu = diag(1./[1,0.35,6]);

    MPC_vars.invTx = diag([3,3,2*pi,4,2,7,30]);
    MPC_vars.invTu = diag([1,0.35,6]);

    MPC_vars.TDu = eye(3);
    MPC_vars.invTDu = eye(3);
    % identity matricies if inputs should not be normalized
    % MPC_vars.Tx = eye(7);
    % MPC_vars.Tu = eye(3);
    % 
    % MPC_vars.invTx = eye(7);
    % MPC_vars.invTu = eye(3);
    % 
    % MPC_vars.TDu = eye(3);
    % MPC_vars.invTDu = eye(3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % state-input bounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % bounds for non-nomalized state-inputs
    % MPC_vars.bounds = [-3,-3,-10,-0.1,-2,-7,   0,  -0.1,-0.35,0  ,  -1 ,-1,-5;
    %                     3, 3, 10,   4, 2, 7,  30,     1, 0.35,5  ,   1 , 1, 5]'; 
    % bounds for nomalized state-inputs (bounds can be changed by changing
    % % normalization)
    MPC_vars.bounds = [-1,-1,-3, 0,-1,-1,   0,  -0.1,-1,0  ,  -1 ,-1,-5;
                        1, 1, 3, 1, 1, 1,   1,     1, 1,1  ,   1 , 1, 5]'; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cost Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MPC_vars.qC = 0.1; % contouring cost
    MPC_vars.qCNmult= 10; % increase of terminal contouring cost
    MPC_vars.qL= 1000; % lag cost
    MPC_vars.qVtheta= 0.02; % theta maximization cost
    MPC_vars.qOmega = 1e-5; % yaw rate regularization cost
    MPC_vars.qOmegaNmult = 10; % yaw rate regularization cost

    MPC_vars.rD= 1e-4; % cost on duty cycle (only used as regularization terms)
    MPC_vars.rDelta= 1e-4; % cost on steering 
    MPC_vars.rVtheta= 1e-4; % cost on virtual velocity

    MPC_vars.rdD= 0.01; % cost on change of duty cycle
    MPC_vars.rdDelta= 1; % cost on change of steering
    MPC_vars.rdVtheta= 0.001; % cost on change of virtual velocity


    MPC_vars.q_eta = 250; % cost on soft constraints (TODO implement soft constraints)

    MPC_vars.costScale = 1; % scaling of the cost for better numerics
    
elseif CarModel == "FullSize"
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MPC settings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % prediction horizon
    MPC_vars.N = 120;
    % sampling time
    MPC_vars.Ts = 0.05;
    % used model (TODO incorparate new models)
    MPC_vars.ModelNo = 2;
    % use bounds on all opt variables (TODO implement selective bounds)
    MPC_vars.fullBound = 1;  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % state-input scaling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % normalization matricies (scale states and inputs to ~ +/- 1
    MPC_vars.Tx = diag(1./[1,1,2*pi,10,10,5,10]);
    MPC_vars.Tu = diag(1./[1,0.5,10]);

    MPC_vars.invTx = diag([1,1,2*pi,10,10,5,10]);
    MPC_vars.invTu = diag([1,0.5,10]);

    MPC_vars.TDu = eye(3);
    MPC_vars.invTDu = eye(3);
    % identity matricies if inputs should not be normalized
    % MPC_vars.Tx = eye(7);
    % MPC_vars.Tu = eye(3);
    % 
    % MPC_vars.invTx = eye(7);
    % MPC_vars.invTu = eye(3);
    % 
    % MPC_vars.TDu = eye(3);
    % MPC_vars.invTDu = eye(3);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % state-input bounds %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % bounds for non-nomalized state-inputs
    % MPC_vars.bounds = [-3,-3,-10,-0.1,-2,-7,   0,  -0.1,-0.35,0  ,  -1 ,-1,-5;
    %                     3, 3, 10,   4, 2, 7,  30,     1, 0.35,5  ,   1 , 1, 5]'; 
    % bounds for nomalized state-inputs (bounds can be changed by changing
    % % normalization)
    MPC_vars.bounds = [-1e4,-1e4,-3, 0.25,-3,-1,   0,    -1,-1, 0  ,  -0.25 ,-0.1,-10;
                        1e4, 1e4, 3,   10, 3, 1, 1e4,     1, 1,10  ,   0.25 , 0.1, 10]'; 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cost Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MPC_vars.qC = 0.01; % contouring cost
    MPC_vars.qCNmult= 10000; % increase of terminal contouring cost
    MPC_vars.qL= 1000; % lag cost
    MPC_vars.qVtheta= 0.5; % theta maximization cost
    MPC_vars.qOmega = 5e0; % yaw rate regularization cost
    MPC_vars.qOmegaNmult = 1; % yaw rate regularization cost

    MPC_vars.rD= 1e-4; % cost on duty cycle (only used as regularization terms)
    MPC_vars.rDelta= 1e-4; % cost on steering 
    MPC_vars.rVtheta= 1e-6; % cost on virtual velocity

    MPC_vars.rdD= 0.1; % cost on change of duty cycle
    MPC_vars.rdDelta= 1; % cost on change of steering
    MPC_vars.rdVtheta= 0.1; % cost on change of virtual velocity

    MPC_vars.q_eta = 250; % cost on soft constraints (TODO implement soft constraints)

    MPC_vars.costScale = 0.01; % scaling of the cost for better numerics
else
    error('invalid model name')
end


end
