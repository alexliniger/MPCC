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

function ModelParams=getModelParams(ModelNo)

if ModelNo==1

    ModelParams.ModelNo=1;
    ModelParams.Scale=1;%scale of the car (1 is a 1:43 scale car)
    
    ModelParams.sx=7; %number of states
    ModelParams.su=3; %number of inputs
    ModelParams.nx=7; %number of states
    ModelParams.nu=3; %number of inputs
    
    ModelParams.stateindex_x=1; %x position
    ModelParams.stateindex_y=2; %y position
    ModelParams.stateindex_phi=3; %orientation
    ModelParams.stateindex_vx=4; %longitudinal velocity
    ModelParams.stateindex_vy=5; %lateral velocity
    ModelParams.stateindex_omega=6; %yaw rate
    ModelParams.stateindex_theta=7; %virtual position

    ModelParams.inputindex_D=1; %duty cycle
    ModelParams.inputindex_delta=2; %steering angle
    ModelParams.inputindex_vtheta=3; %virtual speed
    
    ModelParams.m = 0.041;
    ModelParams.Iz = 27.8e-6;
    ModelParams.lf = 0.029;
    ModelParams.lr = 0.033;

    ModelParams.Cm1=0.287;
    ModelParams.Cm2=0.0545;
    ModelParams.Cr0=0.0518;
    ModelParams.Cr2=0.00035;

    ModelParams.Br = 3.3852;
    ModelParams.Cr = 1.2691;
    ModelParams.Dr = 0.1737;

    ModelParams.Bf = 2.579;
    ModelParams.Cf = 1.2;
    ModelParams.Df = 0.192;
    
    ModelParams.L = 0.12;
    ModelParams.W = 0.06;
    
elseif ModelNo==2
    ModelParams.ModelNo=2;
    ModelParams.Scale=43;%scale of the car (1 is a 1:43 scale car)
    
    ModelParams.sx=7; %number of states
    ModelParams.su=3; %number of inputs
    ModelParams.nx=7; %number of states
    ModelParams.nu=3; %number of inputs
    
    ModelParams.stateindex_x=1; %x position
    ModelParams.stateindex_y=2; %y position
    ModelParams.stateindex_phi=3; %orientation
    ModelParams.stateindex_vx=4; %longitudinal velocity
    ModelParams.stateindex_vy=5; %lateral velocity
    ModelParams.stateindex_omega=6; %yaw rate
    ModelParams.stateindex_theta=7; %virtual position

    ModelParams.inputindex_D=1; %duty cycle
    ModelParams.inputindex_delta=2; %steering angle
    ModelParams.inputindex_vtheta=3; %virtual speed
    
    ModelParams.m = 1573;
    ModelParams.Iz = 2873;
    ModelParams.lf = 1.35;
    ModelParams.lr = 1.35;
    
    ModelParams.Wight_f = ModelParams.lr/(ModelParams.lf+ModelParams.lr);
    ModelParams.Wight_r = ModelParams.lf/(ModelParams.lf+ModelParams.lr);
    
    ModelParams.Cm1=17303;
    ModelParams.Cm2=175;
    ModelParams.Cr0=120;
    ModelParams.Cr2=0.5*1.225*0.35*2.5;%0.5*rho*cd*A

    ModelParams.Br = 13;
    ModelParams.Cr = 2;
    ModelParams.Dr = ModelParams.Wight_f*ModelParams.m*9.81*1.2;

    ModelParams.Bf = 13;
    ModelParams.Cf = 2;
    ModelParams.Df = ModelParams.Wight_r*ModelParams.m*9.81*1.2;
    
    ModelParams.L = 5;
    ModelParams.W = 2.5;
else
    error('ModelNo invalid');
end

end
