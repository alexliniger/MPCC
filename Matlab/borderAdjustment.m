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

function [track,track2] = borderAdjustment(track2,ModelParams,safetyScaling)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% shrink track on both sides with 'WidthCar' %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
WidthCar = 0.5*ModelParams.W*safetyScaling;
%scale track to car size
track2.center = track2.center*ModelParams.Scale;
track2.inner = track2.inner*ModelParams.Scale;
track2.outer = track2.outer*ModelParams.Scale;

% compute width of track (assumption track as a uniform width
widthTrack = norm([track2.inner(1,1)-track2.outer(1,1),track2.inner(2,1)-track2.outer(2,1)]);

for i = 1:length(track2.outer)
    x1 = track2.outer(1,i);
    y1 = track2.outer(2,i);
    x2 = track2.inner(1,i);
    y2 = track2.inner(2,i);
    % vector connecting right and left boundary
    numer=(x2-x1);
    denom=(y1-y2);
    
    % shrinking ratio
    c =  WidthCar/widthTrack;
    d = -WidthCar/widthTrack;
    
    % shrink track
    track.outer(1,i) = x1 + c*numer;
    track.inner(1,i) = x2 - c*numer;
    track.outer(2,i) = y1 + d*denom;
    track.inner(2,i) = y2 - d*denom;    
end
track.center = track2.center;

end