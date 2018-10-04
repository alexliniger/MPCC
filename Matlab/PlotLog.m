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

function [  ] = PlotLog( X,U,Y,track,track2,simN,Ts)
    
    figure(1);
    plot(track.outer(1,:),track.outer(2,:),'r')
    hold on
    plot(track.inner(1,:),track.inner(2,:),'r')
    plot(track2.outer(1,:),track2.outer(2,:),'k')
    plot(track2.inner(1,:),track2.inner(2,:),'k')
    plot(X(1,:),X(2,:),'b')
    if ~isempty(Y)
        for i=1:size(Y,2)
            carBox(Y(:,i),0.025,0.05)
        end
    end
    xlabel('X [m]')
    ylabel('Y [m]')
    axis equal
    hold off

    figure(2)
    subplot(3,1,1)
    plot([0:simN-1]*Ts,U(1,:))
    xlabel('time [s]')
    ylabel('D [-]')
    subplot(3,1,2)
    plot([0:simN-1]*Ts,U(2,:))
    xlabel('time [s]')
    ylabel('\delta [rad]')
    subplot(3,1,3)
    plot([0:simN-1]*Ts,U(3,:))
    xlabel('time [s]')
    ylabel('v_{\Theta} [m/s]')

    figure(3)
    subplot(3,1,1)
    plot([0:simN-1]*Ts,X(4,:))
    xlabel('time [s]')
    ylabel('\phi [rad]')
    subplot(3,1,2)
    plot([0:simN-1]*Ts,X(6,:))
    xlabel('time [s]')
    ylabel('\omega [rad/s]')
    subplot(3,1,3)
    plot([0:simN-1]*Ts,180/pi*atan2(X(5,:),X(4,:)))
    xlabel('time [s]')
    ylabel('\beta [deg]')

    pause(0.001)
end

