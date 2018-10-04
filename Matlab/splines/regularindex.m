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

function ind=regularindex(ppx,ppy,number)
umax=ppx.breaks(end);
l=splinelength(ppx, ppy,ppx.breaks(1:(end-1)),ppx.breaks(2:end));
cl=cumsum(l)/sum(l);
ind=zeros(number,1);
ind(1:floor(number*cl(1)))=linspace(ppx.breaks(1),ppx.breaks(2),floor(number*cl(1)));
for i=2:(length(l))
    ind(floor(number*cl(i-1)):floor(number*cl(i)))=linspace(ppx.breaks(i),ppx.breaks(i+1),floor(number*cl(i))-floor(number*cl(i-1))+1);
end

end