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

function [ppx ppy err]=splineInterp(X,Y,Tsample,circular)
%[ppx ppy err]=splineInterp(X,Y,Tsample)
%Compute a 2D spline using X and Y. Optionally compute the error made by sampling (zero if Ts=1)
%Tsample has to be an integer 1<=Tsample<=legnth(X), by default Tsample=1
%if circular is set to 'y', the function add a point to close the
%trajectory and ensure derivability at closure.
%X and Y and 1D data, length(X)=length(Y)
%err represent the squared error due to sampling
capprox=50;%number of part to evaluate the spline when calculating the error
if(nargin<3)
    Tsample=1;
end
if(nargin<4)
    circular='n';
end

t=0:(ceil(length(X)/Tsample)-1);
s=1:Tsample:length(X);
if(length(X)/Tsample~=floor(length(X)/Tsample))
    s=[s length(X)];
    t = [t t(end)+1];
end

if(strcmp(circular,'y')==1)
    t=[t t(end)+1];%add a point to close the trajectory;
    %this is not perfect because the value of the derivative should be
    %included in the unknowns of linear solver used for spline calculation.
    theta=atan2(Y(1)-Y(end-1),X(1)-X(end-1));
    r=0;
    Xs=[X(s) X(s(1))];
    Ys=[Y(s) Y(s(1))];
else
    Xs=X(s);
    Ys=Y(s);
end
% figure(3);plot(Xs,Ys,'xr-');
ppx=spline(t,Xs);
ppy=spline(t,Ys);

if(nargout>2)%calculate the error
    err=0;
    if(Tsample==1)
        return
    end
    for i=t(1:(end-1))%for each piece
        x=ppval(ppx,linspace(i,i+1,capprox));
        y=ppval(ppy,linspace(i,i+1,capprox));
        for j=1:Tsample%for each point
            if(i*floor(length(X)/Tsample)+j>length(X))
                break;
            end
            err= err + min((x-X(i*floor(length(X)/Tsample)+j)).^2+(y-Y(i*floor(length(Y)/Tsample)+j)).^2);
        end
    end
    err = err / length(X);
end

end
