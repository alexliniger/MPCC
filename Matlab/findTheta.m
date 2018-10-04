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

function [ theta, closestIdx] = findTheta(currentPosition,TrackCenter,traj_breaks,trackWidth,last_closestIdx)
%FINDTHETA returns theta, the index of the point on the centerline which is the
%closest to the current position

%the algorithm first performs a local search around a guess where the
%closest point on the center line is if the error is too big (more than
%half the track width) a global search is performed.

%then the projection is computed using the inner product and assuming the
%track is picewise affine

posX = currentPosition(1);
posY = currentPosition(2);
TrackX = TrackCenter(1,:);
TrackY = TrackCenter(2,:);
trackCenterline = TrackCenter;
N_Track = length(TrackX);

% length of local search region
N_searchRegion = 26;
searchRegion = zeros(1,N_searchRegion);
% determin index which are investigated
k = 1;
for i = last_closestIdx-5:last_closestIdx+20
    searchRegion(k) = i;
    k = k+1;
end
% wrapping if necessary
if last_closestIdx > N_Track - 20 || last_closestIdx <= 5
    i = find(searchRegion < 1);
    k = find(searchRegion > N_Track);
    searchRegion(i) = N_Track+searchRegion(i);
    searchRegion(k) = searchRegion(k)-N_Track;
end

% compute squared distance to the current location
TrackXsmall = TrackX(searchRegion);
TrackYsmall = TrackY(searchRegion);
distanceX = TrackXsmall - posX*ones(1,N_searchRegion);
distanceY = TrackYsmall - posY*ones(1,N_searchRegion);

squared_distance_array   = distanceX.^2 + distanceY.^2;
% find closest point
[e , minIndex] = min(squared_distance_array);
minIndex = searchRegion(minIndex);

% if distance is to large perform global search with respect to track width
if sqrt(e) > (trackWidth*1.25)/2
    distanceX2 = TrackX - posX*ones(1,N_Track);
    distanceY2 = TrackY - posY*ones(1,N_Track);

    squared_distance_array2   = distanceX2.^2 + distanceY2.^2;

    [e , minIndex] = min(squared_distance_array2);
end

% circular edge conditions
if minIndex == 1;
    nextIdx = 2;
    prevIdx = N_Track;
elseif minIndex == N_Track;
    nextIdx = 1;
    prevIdx = N_Track-1;
else
    nextIdx = minIndex + 1;
    prevIdx = minIndex - 1;
end


% compute theta based on inner product projection
closestIdx = minIndex;
cosinus = dot([posX;posY] - trackCenterline(:,minIndex), trackCenterline(:,prevIdx) - trackCenterline(:,minIndex));
if (cosinus > 0)
    minIndex2 = prevIdx;
else
    minIndex2 = minIndex;
    minIndex = nextIdx;
end


if e ~= 0
    cosinus = dot([posX;posY] - trackCenterline(:,minIndex2), trackCenterline(:,minIndex) - trackCenterline(:,minIndex2))/...
                 (norm([posX;posY] - trackCenterline(:,minIndex2))*norm(trackCenterline(:,minIndex) - trackCenterline(:,minIndex2)));
else
    cosinus = 0;
end
theta = traj_breaks(1,minIndex2);
theta = theta + cosinus*norm([posX;posY] - trackCenterline(:,minIndex2),2);
return