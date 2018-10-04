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

function [ppcx,ppcy]=computeCenter(ppx,ppy,Xtr,Ytr,Tsample)
% Compute the spline approximation of the center of the track.
%
% Works by subdividing every two intervals, splining that interval, and
% defining a new spline point as the point closest to the splined
% subinterval.
%
% Modified by Samuel Zhao Oct.18,2012 so borders is same length as traj
% This function assumes circular end-condition for track.


refine=50;%factor of refinement to compute correspondance
nbreaks=numel(ppx.breaks);
Xtrs=[Xtr(1:end) Xtr(end)]; % wrap-around
Ytrs=[Ytr(1:end) Ytr(end)]; % wrap-around
t=zeros(nbreaks,1);
Xj=ppval(ppx,ppx.breaks);
Yj=ppval(ppy,ppx.breaks);

t(1)=ppx.breaks(1); %enforce start and end points are same
t(end)=ppx.breaks(end);

for i=2:nbreaks-1

    cX=Xtrs(i);
    cY=Ytrs(i);

    [~,b]=min((Xj-cX).^2+(Yj-cY).^2);%find closest point
    
    % Pnext is spline parameter of next track coordinate
    % Pprev is spline parameter of previous track coordinate
    if(b==nbreaks)
        Pnext=ppx.breaks(b);
    else
        Pnext=ppx.breaks(b+1);
    end
    if(b==1)
        Pprev=ppx.breaks(1);
    else
        Pprev=ppx.breaks(b-1);
    end
    Xref=ppval(ppx,linspace(Pprev,Pnext,refine));
    Yref=ppval(ppy,linspace(Pprev,Pnext,refine));
    [~,d]=min((Xref-cX).^2+(Yref-cY).^2);%find closest point
%    bla(i)=Pnext;
    index=Pprev+d*(Pnext-Pprev)/refine;
    t(i)=index;

end

ppcx=spline(t,Xtrs);
ppcy=spline(t,Ytrs);
end
