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

function config=getMpcConfig() 
    config.sx=11; %number of states
    config.su=4; %number of inputs
    config.nx=11; %number of states
    config.nu=4; %number of inputs
    
    config.stateindex_x=1; %x position
    config.stateindex_y=2; %y position
    config.stateindex_yaw=3; %orientation
    config.stateindex_vx=4; %longitudinal velocity
    config.stateindex_vy=5; %lateral velocity
    config.stateindex_r=6; %yaw rate
    config.stateindex_s=7; %virtual position
    config.stateindex_throttle=8; %throttle
    config.stateindex_steeringAngle=9; %steering angle
    config.stateindex_brakes=10; %brakes
    config.stateindex_vs=11; %virtual speed

    config.inputindex_dThrottle=1; %throttle rate of change
    config.inputindex_dSteeringAngle=2; %steering angle rate of change
    config.inputindex_dBrakes=3; %brakes rate of change
    config.inputindex_dVtheta=4; %virtual speed rate of change
end
