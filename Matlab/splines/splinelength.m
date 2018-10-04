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

function l=splinelength(ppx, ppy,umin,umax)
%l=splinelength(ppx, ppy,umin,umax)
%compute the length of a 2D spline between parameter umin and umax
%ppx and ppy are spline parameter. See spline function for more details
%umin and umax are vectors of the same length
%and for all i, umin(i)<=umax(i)
%checked
assert(length(umin)==length(umax));

dppx=getSplineDerivatives(ppx);
dppy=getSplineDerivatives(ppy);

h = @(u) sqrt(ppval(dppx,u).^2+ppval(dppy,u).^2);
l=zeros(length(umin),1);
for i=1:length(umin)
    l(i)=quad(h,umin(i),umax(i));
end

end
