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

function xp=SimTimeStep(x,u,Ts,ModelParams)
%x state
%u input
%Ts sampling time
x0=x;

[~,inivt]=ode45(@(t,x)fx_bicycle(t,x,u,ModelParams),[0 Ts],x0);
xp=inivt(end,:);
return

function xdot=fx_bicycle(t,x,u,ModelParams)

    Cm1=ModelParams.Cm1;
    Cm2=ModelParams.Cm2;
    Cr0=ModelParams.Cr0;
    Cr2=ModelParams.Cr2;
    
            
    B_r = ModelParams.Br;
    C_r = ModelParams.Cr;
    D_r = ModelParams.Dr;
    B_f = ModelParams.Bf;
    C_f = ModelParams.Cf;
    D_f = ModelParams.Df;
    
    m = ModelParams.m;
    Iz = ModelParams.Iz;
    l_f = ModelParams.lf;
    l_r = ModelParams.lr;

    
    phi   =x(3);
    v_x     =x(4);
    v_y     =x(5);
    omega   =x(6);
    D     =u(1);
    delta =u(2);
    
    
    alpha_f = -atan2(l_f*omega + v_y,abs(v_x))+delta;
    alpha_r =  atan2(l_r*omega - v_y,abs(v_x));

    F_fy = D_f*sin(C_f*atan(B_f*alpha_f));
    F_ry = D_r*sin(C_r*atan(B_r*alpha_r));

    F_rx = (Cm1*D-Cm2*D*v_x-Cr0-Cr2*v_x^2);

    xdot=[v_x*cos(phi) - v_y*sin(phi);
       v_y*cos(phi) + v_x*sin(phi);
       omega;
       1/m*(F_rx - F_fy*sin(delta) + m*v_y*omega);
       1/m*(F_ry + F_fy*cos(delta) - m*v_x*omega);
       1/Iz*(F_fy*l_f*cos(delta)- F_ry*l_r);
       u(3)];

    
return