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

function [nppx, nppy, err]=normalizedSplineInterp(X,Y,Tsample,circular)
%[nppx, nppy, err]=normalizedSplineInterp(X,Y,Tsample)
%return a quasi-arc-parametrized spline
%Tsample is optional and default value is 1
%if circular is set to 'y', the function add a point to close the
%trajectory and ensure derivability at closure.
%checked for cubic splines
if(nargin<3)
    Tsample=1;
end
if(nargin<4)
    circular='n';
end
if(nargout>2)
    [nppx nppy err]=splineInterp(X,Y,Tsample,circular);
else
    [nppx nppy]=splineInterp(X,Y,Tsample,circular);
end
l=splinelength(nppx, nppy,nppx.breaks(1:(end-1)),nppx.breaks(2:end));
cl=cumsum(l);

%vandermonde matrix
h=@(x) [x.^3 x.^2 x ones(numel(x),1)];
div=h(1./l);

nppx.breaks = [0 cl'];
nppy.breaks = [0 cl'];
nppx.coefs = nppx.coefs .* div;
nppy.coefs = nppy.coefs .* div;
end