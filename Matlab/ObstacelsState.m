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

function Y = ObstacelsState(track,traj,trackWidth,n_cars)

if n_cars == 1
    Y = [];
elseif n_cars == 5
    startIdx_Op1 = 620; 
    vx0_Op1 = 2; 
    % find theta that coresponds to the 10th point on the centerline
    [theta_Op1, ~] = findTheta([track.center(1,startIdx_Op1),track.center(2,startIdx_Op1)],track.center,traj.ppx.breaks,trackWidth,startIdx_Op1);

    x0_Op1 = [track.center(1,startIdx_Op1),track.center(2,startIdx_Op1)+trackWidth/2,... % point on centerline
          atan2(ppval(traj.dppy,theta_Op1),ppval(traj.dppx,theta_Op1)),... % aligned with centerline
          vx0_Op1 ,0,0,theta_Op1]'; %driving straight with vx0, and correct theta progress

    startIdx_Op2 = 625; 
    vx0_Op2 = 2; 
    % find theta that coresponds to the 10th point on the centerline
    [theta_Op2, ~] = findTheta([track.center(1,startIdx_Op2),track.center(2,startIdx_Op2)],track.center,traj.ppx.breaks,trackWidth,startIdx_Op2);

    x0_Op2 = [track.center(1,startIdx_Op2),track.center(2,startIdx_Op2)+trackWidth/4,... % point on centerline
          atan2(ppval(traj.dppy,theta_Op2),ppval(traj.dppx,theta_Op2)),... % aligned with centerline
          vx0_Op2 ,0,0,theta_Op2]'; %driving straight with vx0, and correct theta progress

    startIdx_Op3 = 630; 
    vx0_Op3 = 2; 
    % find theta that coresponds to the 10th point on the centerline
    [theta_Op3, ~] = findTheta([track.center(1,startIdx_Op3),track.center(2,startIdx_Op3)],track.center,traj.ppx.breaks,trackWidth,startIdx_Op3);

    x0_Op3 = [track.center(1,startIdx_Op3),track.center(2,startIdx_Op3),... % point on centerline
          atan2(ppval(traj.dppy,theta_Op3),ppval(traj.dppx,theta_Op3)),... % aligned with centerline
          vx0_Op3 ,0,0,theta_Op3]'; %driving straight with vx0, and correct theta progress

    startIdx_Op4 = 645; 
    vx0_Op4 = 2; 
    % find theta that coresponds to the 10th point on the centerline
    [theta_Op4, ~] = findTheta([track.center(1,startIdx_Op4),track.center(2,startIdx_Op4)],track.center,traj.ppx.breaks,trackWidth,startIdx_Op4);

    x0_Op4 = [track.center(1,startIdx_Op4),track.center(2,startIdx_Op4)-trackWidth/3,... % point on centerline
          atan2(ppval(traj.dppy,theta_Op4),ppval(traj.dppx,theta_Op4)),... % aligned with centerline
          vx0_Op4 ,0,0,theta_Op4]'; %driving straight with vx0, and correct theta progress
    
    Y = [x0_Op1,x0_Op2,x0_Op3,x0_Op4];
else
    error('only 0 or 4 obstacles are pre programmed settings, feel free to change the obstacle constellation')
end