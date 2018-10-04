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

function [  ] = PlotPrediction( x,u,b,Y,track,track2,traj,MPC_vars,ModelParams )
    
    N = MPC_vars.N;
    Ts = MPC_vars.Ts;
    
    tl = traj.ppx.breaks(end);
    
    figure(1);
    plot(track.outer(1,:),track.outer(2,:),'r')
    hold on
    plot(track.inner(1,:),track.inner(2,:),'r')
    plot(track2.outer(1,:),track2.outer(2,:),'k')
    plot(track2.inner(1,:),track2.inner(2,:),'k')
    plot(b(:,1),b(:,2),'g.')
    plot(b(:,3),b(:,4),'g.')
    plot(ppval(traj.ppx,mod(x(7,:),tl)),ppval(traj.ppy,mod(x(7,:),tl)),':k')
    plot(x(1,:),x(2,:),'-b')
    carBox(x(:,1),ModelParams.W/2,ModelParams.L/2)
    if ~isempty(Y)
        for i=1:size(Y,2)
            carBox(Y(:,i),ModelParams.W/2,ModelParams.L/2)
        end
    end
    xlabel('X [m]')
    ylabel('Y [m]')
    axis equal
    hold off

    figure(2)
    subplot(3,1,1)
    plot([0:N-1]*Ts,u(1,:))
    xlabel('time [s]')
    ylabel('D [-]')
    subplot(3,1,2)
    plot([0:N-1]*Ts,u(2,:))
    xlabel('time [s]')
    ylabel('\delta [rad]')
    subplot(3,1,3)
    plot([0:N-1]*Ts,u(3,:))
    xlabel('time [s]')
    ylabel('v_{\Theta} [m/s]')

    figure(3)
    subplot(3,1,1)
    plot([0:N]*Ts,x(3,:))
    xlabel('time [s]')
    ylabel('\phi [rad]')
    subplot(3,1,2)
    plot([0:N]*Ts,x(6,:))
    xlabel('time [s]')
    ylabel('\omega [rad/s]')
    subplot(3,1,3)
    plot([0:N]*Ts,x(7,:))
    xlabel('time [s]')
    ylabel('\Theta [m]')

    pause(0.001)

end

