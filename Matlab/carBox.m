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

function [] = carBox(x0,w,l)
        car1 = x0(1:2) + [cos(x0(3))*l;sin(x0(3))*l] + [sin(x0(3))*w;-cos(x0(3))*w];
        car2 = x0(1:2) + [cos(x0(3))*l;sin(x0(3))*l] - [sin(x0(3))*w;-cos(x0(3))*w];
        car3 = x0(1:2) - [cos(x0(3))*l;sin(x0(3))*l] + [sin(x0(3))*w;-cos(x0(3))*w];
        car4 = x0(1:2) - [cos(x0(3))*l;sin(x0(3))*l] - [sin(x0(3))*w;-cos(x0(3))*w];

        plot([car1(1),car2(1),car4(1),car3(1),car1(1)],[car1(2),car2(2),car4(2),car3(2),car1(2)],'k','LineWidth',1)
end
