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

function [new_b_array, last_border_array]= getNewBorders(traj,borders,track_center,X_all,n_cars,caridx,MPC_vars,ModelParams)
% this function will replace the carinsight function
%
% Inputs:
%   traj        - pathinfo (gotten from splinify)
%   borders     - trackinfo (gotten from splinify)
%   track_center- track.center (where track is what went through splinify)
%   X_all       - state vector of all cars
%   n_cars      - total number of cars
%   carindex    - index of controlled car
%   MPC_vars    - MPC parameters
%   ModelParams - Model parameters
%
% Outputs:
%   new_b_array - [N]x[4]x[n_cars] matrix containing the points to define
%                 corridors for each car
%  last_border_array - [ncars]x[3] matrix containing the points to define
%                 corridors for each car
%
% Pseudo code:
%
% for each car (caridx)
%   for all other cars (othercaridx)
%       if othercaridx is "close enough" to caridx
%           populate grid 
%       end
%   end
% end
% 
% Note: X is length (N+1)*nx, while grid has N rows. The i'th row of grid
% corresponds to the (i+1)'th section of X. That is, x0 does NOT have a
% corresponding line in grid.

% design parameters
DesignParameters.N = MPC_vars.N;
DesignParameters.Ts = MPC_vars.Ts;
DesignParameters.grid_width = 12; % track width is split into grid_width+1 equal sections

%   with margin
DesignParameters.carwidth =  ModelParams.W*2; %[m]
DesignParameters.carlength = ModelParams.L*2; %[m] 

lookahead = X_all(MPC_vars.N*ModelParams.nx+ModelParams.stateindex_theta,1) - X_all(ModelParams.stateindex_theta,1);
DesignParameters.thetadifferencethreshold=lookahead*1.25; %[m] (used by CarIsClose() )
DesignParameters.Cost_dead_end = 5000;
DesignParameters.Cost_maxanglechange = 2000; % cost when maxanglechange is exceeded. Choose a cost that is greater than Cost_dead_end. used in getCost_AngleChange()%
DesignParameters.Cost_Weights = [0.5 3 0.5]; % [deviation length anglechange]
DesignParameters.time_step_cost_discount_factor = 0.95;

% misc. constants
Constants.UNOCCUPIED=0; % fill value of grid_isOccupiedBy if unoccupied (used by getOptimalPathMatrix() and Find_Nearest_Borders() )
Constants.DEAD_END_PATH = -2; % fill value of Optm_Path_Mtx if no valid paths (used by getOptimalPathMatrix() and Find_Nearest_Borders() )
Constants.LAST_BORDER_XMAX = 10.0; % [m] largest value for x, used as default last_border

nx = ModelParams.nx;
N = MPC_vars.N;

% %%%%%%%%%%%%% end of control variables %%%%%%%%%%%%%%%%%%%%

    grid_width = DesignParameters.grid_width;

    % misc. defines
    tracklength = traj.ppx.breaks(end);

    % declare variable
    % new_b_array=zeros(N,4); %[N]x[4]x[n_cars] matrix 
    last_border_array=zeros(3,n_cars); %[3]x[n_cars] matrix 


    % grid coordinates for current car
    X = X_all(:,caridx); % full state vector of current car
    [grid_x, grid_y, new_b_default] = getGridCoordinates(traj,borders,track_center,X,DesignParameters,ModelParams);

    grid_isOccupiedBy = zeros(N,grid_width); %initialize
    % for c, type cast this to integers
    % see documents: mxGetData, mxClassID and cast
    grid_isOccupiedBy = cast(grid_isOccupiedBy, 'int32');

    %Populate grid_isOccupiedBy
    for othercaridx=1:n_cars
        if othercaridx==caridx,continue;end
        
        X_oth = X_all(:,othercaridx);

        %Check to see if car is even remotely close
        if CarIsClose(X,X_oth,nx,tracklength,DesignParameters)
            %fill relevant values of grid_isOccupiedBy, if any
            grid_isOccupiedBy = getGridOccupancy(grid_isOccupiedBy,grid_x,grid_y,X_oth,othercaridx,...
                        DesignParameters,ModelParams);
        end
    end


    OptmPath = 0;

    % new_b_array{caridx} is a [N]by[4] matrix. Each row contains two
    % (x,y) coordinates marking the corridor
    if max(max(abs(grid_isOccupiedBy)))==0 %in C -1, first car is car0
        % if all elements of grid_isOccupiedBy are zeros (ie. no cars nearby),
        % then use default borders
        new_b_array = new_b_default;
        last_border_array = [1;0;10]; % x<=10
    else
%         try
        [new_b_array last_border_array OptmPath] = getNewBorders_helper(X_all,X,grid_x,grid_y,new_b_default,grid_isOccupiedBy,DesignParameters,ModelParams,Constants);

        %%PLOTS 
        boxSize = [DesignParameters.carlength, DesignParameters.carwidth];
        % plotDPdata(X_all, caridx, boxSize, grid_x, grid_y, grid_isOccupiedBy, OptmPath);

    end
    


% %%PLOTS 
% boxSize = [DesignParameters.carlength, DesignParameters.carwidth];
% plotDPdata(X_all, caridx, boxSize, grid_x, grid_y, grid_isOccupiedBy, OptmPath);
end



function [grid_x grid_y new_b_default] = getGridCoordinates(traj,borders,track_center,X,DesignParameters,ModelParams)
% assumes x,y are the 1 and 2 states



nx=ModelParams.nx;
N=DesignParameters.N;
grid_width=DesignParameters.grid_width;

tracklength = traj.ppx.breaks(end);

% find physical thetas for each future position in X
theta_phys=zeros(N,1);
for i=1:N
    xy_posn=X(i*nx+1:i*nx+2);% assume x,y are the 1 and 2 states
    theta_virt = X(nx*i+nx);

    % approximate theta_phys as theta_virt: should hold if lag error weighting is high enough
    theta_virt=mod(theta_virt,tracklength);
    theta_phys(i)=theta_virt; 

end

TrackLeftx=ppval(borders.pplx,theta_phys);
TrackLefty=ppval(borders.pply,theta_phys);
TrackRightx=ppval(borders.pprx,theta_phys);
TrackRighty=ppval(borders.ppry,theta_phys);
Dx = TrackRightx - TrackLeftx;
Dy = TrackRighty - TrackLefty;

% border if no car is near
new_b_default = [TrackLeftx TrackLefty TrackRightx TrackRighty];

grid_x = zeros(N,grid_width);
grid_y = zeros(N,grid_width);
for i=1:N
    for j=1:grid_width
        % equally spaced grid
        grid_x(i,j) = TrackLeftx(i) + j/(grid_width+1)*Dx(i);
        grid_y(i,j) = TrackLefty(i) + j/(grid_width+1)*Dy(i);
    end
end    

end

function bool = CarIsClose(X_curr,X_oth,nx,tracklength,DesignParameters) 
% This implements a quick and cheap method to check if two cars are close.
% Simply compare their theta values at time 0. If they are within a
% threshold of each other, returns true. Otherwise, returns false. Note
% that wrap-around is implemented with respect to tracklength.
%
% Note: assumed theta index is nx

threshold=DesignParameters.thetadifferencethreshold;


bool=false;

if               abs(X_curr(nx) - X_oth(nx)) < threshold || ...
   tracklength - abs(X_curr(nx) - X_oth(nx)) < threshold
    % 2nd line of above if-condition checks for track wrap-around
    bool=true;
    return;
end

end

function grid_isOccupiedBy = getGridOccupancy(grid_isOccupiedBy,grid_x,grid_y,X_oth,othercaridx,...
    DesignParameters,ModelParams)
% check if (grid_x,grid_y)(i,j) is inside polygon of {othercaridx}'th car at
% corresponding time.
%
% if grid(i,j) is occupied by the car, grid_isOccupiedBy(i,j)
% is assigned the index of the occupying car (othercaridx). Otherwise it is
% untouched.

carwidth=DesignParameters.carwidth;
carlength=DesignParameters.carlength;
stateindex_x=ModelParams.stateindex_x;
stateindex_y=ModelParams.stateindex_y;
stateindex_phi=ModelParams.stateindex_phi;

nx=ModelParams.nx;
N=DesignParameters.N;


for i=1:N
    % position of car
    posn = [X_oth(i*nx+stateindex_x) X_oth(i*nx+stateindex_y)];
    
    % orientation of car (zero is positive x axis)
    phi = X_oth(i*nx+stateindex_phi);
    cosphi=cos(phi);
    sinphi=sin(phi);
    
    % tangent and normal vectors of bounding box (normal vector is 90 deg
    % clockwise of tangent vector)
    t = [cosphi,sinphi]*carlength/2;
    n = [sinphi,-cosphi]*carwidth/2;
    
    % corners of bounding box of the car with respect to the point of
    % rotation used in the model
    FR = posn + t + n;
    FL = posn + t - n;
    BR = posn - t + n;
    BL = posn - t - n;
    
    % check gridpoints of i'th row is inside bounding box of car
    vertices = [FR;BR;BL;FL]; % [number of corners]by[2] array
    points = [grid_x(i,:);grid_y(i,:)]'; % [number of points to test]by[2] array

    in = inpoly(points,vertices);
    
    grid_isOccupiedBy(i,in) = othercaridx;

end

end

function [new_b last_border OptmPath] = getNewBorders_helper(X_all,X,grid_x,grid_y,new_b_default,grid_isOccupiedBy,DesignParameters,ModelParams,Constants)
% to do


priorPath = Lock2Grid(X,grid_x,grid_y,DesignParameters,ModelParams);

% get [N]by[grid_width]by[grid_width] matrix of optimal policies
[Optm_Path_Mtx Optm_Cost_Mtx]=getOptimalPathMatrix(X,grid_x,grid_y,grid_isOccupiedBy,priorPath,DesignParameters,ModelParams,Constants);
% Optm_Cost_Mtx
% trace optimal policy matrix to get optimal path
OptmPath=TraceOptimalPath(Optm_Path_Mtx,DesignParameters,Constants);

% find borders of optimal path
[new_b last_border] = Find_Nearest_Borders(X_all,grid_x,grid_y,new_b_default,grid_isOccupiedBy,OptmPath,DesignParameters,ModelParams,Constants);

end

function priorPath = Lock2Grid(X,grid_x,grid_y,DesignParameters,ModelParams)
% priorPath - [N]by[1] (dimensionless)

nx = ModelParams.nx;
N = DesignParameters.N;
grid_width = DesignParameters.grid_width;
stateindex_x=ModelParams.stateindex_x;
stateindex_y=ModelParams.stateindex_y;

priorPath = zeros(N,1);
% for c, type cast this to integers
priorPath = cast(priorPath, 'int32');
    
for i=1:N
    min=100; % choose some big number that distance^2 can never reach
    ind=1;
    for j=1:grid_width
        dx = grid_x(i,j) - X(i*nx+stateindex_x);
        dy = grid_y(i,j) - X(i*nx+stateindex_y);
        d2=dx*dx+dy*dy;
        if d2<min
            min=d2;
            ind=j;
        end
    end
    priorPath(i) = ind;
end

end

function [Optm_Path_Mtx Optm_Cost_Mtx]=getOptimalPathMatrix(X,grid_x,grid_y,grid_isOccupiedBy,priorPath,DesignParameters,ModelParams,Constants)
% for i>1: jnxt=Optm_Path_Mtx(i+1,j,j1) is the optimal next position
%          (i+1,jnxt) if you are at (i,j) and you came from (i-1,j1)
%          (goes up to i=N-1)
% for i=1: Optm_Path_Mtx(2,j,1) is the optimal next position if you are at
%          (1,j)
% for i=0: Optm_Path_Mtx(1,1,1) is the optimal next position from the
%          starting position
%
% Optm_Cost_Mtx(i+1,j,j1) is the optimal cost if you are at (i,j) and came
% from (i-1,j1) and choose to go to (i+1,jnxt)
%
% note: if (i+1,jnxt) is invalid for all jnxt, then Optm_Path_Mtx(i+1,j,j1)
% is assigned value of DEAD_END_PATH. Note that the position (i,j) itself
% is a valid position.

run_getOptimalPathMatrix_mex=0;
if run_getOptimalPathMatrix_mex
[Optm_Path_Mtx_mex Optm_Cost_Mtx_mex] = getOptimalPathMatrix_mex(X,grid_x,grid_y,grid_isOccupiedBy,priorPath);
end

N=DesignParameters.N;
grid_width=DesignParameters.grid_width;
Cost_dead_end = DesignParameters.Cost_dead_end;
time_step_cost_discount_factor = DesignParameters.time_step_cost_discount_factor;

% nx=ModelParams.nx;
stateindex_x=ModelParams.stateindex_x;
stateindex_y=ModelParams.stateindex_y;
% stateindex_phi=ModelParams.stateindex_phi;
% stateindex_v=ModelParams.stateindex_v;

UNOCCUPIED=Constants.UNOCCUPIED;
DEAD_END_PATH = Constants.DEAD_END_PATH;



Optm_Path_Mtx=zeros(N,grid_width,grid_width);
Optm_Path_Mtx = cast(Optm_Path_Mtx, 'int32');


% figure(4); title('getNewBorders.m\getOptimalPathMatrix()\grid_isOccupiedBy');
% spy(grid_isOccupiedBy)
% asdf=5;

for i=N-1:-1:1
    for j1=1:grid_width
        % no need to look for a valid path if previous position (i-1,j1) is
        % invalid. Also no need to check for i=1, since i=0 will always
        % have a valid j
        if i>1 && grid_isOccupiedBy(i-1,j1)~=UNOCCUPIED
            continue;
        end
        
        % change to true if *any* valid jnxt is found for at least 1 j
        valid_path_is_found = false; 
        
        for j=1:grid_width
            % no need to look for a valid path if current position (i,j) is
            % invalid
            if grid_isOccupiedBy(i,j)~=UNOCCUPIED
                continue;
            end

            
            % change to true if a valid jnxt is found for *current* j
            valid_jnxt_is_found = false;
        

            for jnxt=1:grid_width
                %if (i<=N-2 && Optm_Path_Mtx(i+2,jnxt,j)==DEAD_END_PATH) || ...
                if    PathIsObstructed(grid_isOccupiedBy,i,j,jnxt,UNOCCUPIED)
                    
                    if jnxt==grid_width && ~valid_jnxt_is_found
                        %uh oh: no valid jnxt found for this j
                        %assign a really big cost to this state
                        Optm_Cost_Mtx(i+1,j,j1) = Cost_dead_end;
                        Optm_Path_Mtx(i+1,j,j1) = DEAD_END_PATH;

                    end
                    
                    % if not yet at last point in row, continue checking.
                    % if at last point, move onto the next j.
                    continue;
                end
                
                % if code reaches here, both j and jnxt are valid

                % calculate costs
                Cost_PathDeviation=getCost_PathDeviation(priorPath,i,jnxt);
                Cost_PathLength=getCost_PathLength(grid_x,grid_y,X,i,j,jnxt,stateindex_x,stateindex_y);
                Cost_AngleChange=getCost_AngleChange(grid_x,grid_y,X,i,j,j1,jnxt,DesignParameters,ModelParams);
                
                % add costs
                if i==N-1
                    Cost_Total = DesignParameters.Cost_Weights*[Cost_PathDeviation;Cost_PathLength;Cost_AngleChange];
                else
                    % also add on (discounted) cost of (i+2,jnxt,j)
                    Cost_Total = DesignParameters.Cost_Weights*[Cost_PathDeviation;Cost_PathLength;Cost_AngleChange] ...
                        + time_step_cost_discount_factor * Optm_Cost_Mtx(i+2,jnxt,j);
                end
                
%                 if(i==12 && j==4 && j1==1)
%                     tmp=[i j j1 jnxt;
%                         Cost_PathDeviation Cost_PathLength Cost_AngleChange Optm_Cost_Mtx(i+2,jnxt,j);
%                         Cost_Total zeros(1,3)];
%                     asdf=5;
%                 end
                
                if ~valid_jnxt_is_found
                    % this is first valid path found, so set values for
                    % later comparison
                    minCost = Cost_Total;
                    jnxt_minCost = jnxt;
                    valid_jnxt_is_found = true;
                else
                    % compare with previous minimum cost
                    if Cost_Total<minCost
                        minCost = Cost_Total;
                        jnxt_minCost=jnxt;
                    end
                end
            end %end of jnxt loop

            
            

            if valid_jnxt_is_found
                % minCost is now the minimum cost for (i,j,j1)
                Optm_Cost_Mtx(i+1,j,j1) = minCost;
                Optm_Path_Mtx(i+1,j,j1) = jnxt_minCost;
                valid_path_is_found = true;
            end
                       
            
        end % end of j loop
        
        if ~valid_path_is_found
            %uh oh: no valid path is found!
            % return an error code and exit
            error('no_valid_path');
        end
        
        if i==1
            % for i=1, does not use j1, so arbitrarily choose j1=1 to
            % store data in. Hence, can break out of the j1 loop.
            break;
        end
    end % end of j1 loop
    
end % end of i loop


% i=0 must be dealt with separately
i=0;
valid_jnxt_is_found=false; %since only the single j, find a valid jnxt IFF find a valid path

for jnxt=1:grid_width
    % check if grid(1,jnxt) is occupied
    if grid_isOccupiedBy(1,jnxt)~=UNOCCUPIED
        if jnxt==grid_width && ~valid_jnxt_is_found
            % uh oh: no valid path for any jnxt!
            % return an error code and exit
            error('no_valid_path');
        end
        % if not yet reached last point in row, continue checking
        continue;
    end
    
    % if code reaches here, valid jnxt has been found
    
    % costs
    Cost_PathDeviation=getCost_PathDeviation(priorPath,i,jnxt);
    Cost_PathLength=getCost_PathLength(grid_x,grid_y,X,i,j,jnxt,stateindex_x,stateindex_y);
    Cost_AngleChange=getCost_AngleChange(grid_x,grid_y,X,i,j,j1,jnxt,DesignParameters,ModelParams); %j1 NOT SET...
    
    % total cost
    Cost_Total = DesignParameters.Cost_Weights*[Cost_PathDeviation;Cost_PathLength;Cost_AngleChange] ...
        + time_step_cost_discount_factor ...
        * Optm_Cost_Mtx(i+2,jnxt,1); % 3rd dim=1 here because j1 for i=1 is always 1
    
    if ~valid_jnxt_is_found
        % first valid path found, so set values for later comparison
        minCost=Cost_Total;
        jnxt_minCost=jnxt;
        valid_jnxt_is_found=true;
    else
        % compare with previous minimum cost
        if Cost_Total<minCost
            minCost = Cost_Total;
            jnxt_minCost=jnxt;
        end
    end
    
end % end of jnxt loop

% minCost is now the minimum cost for Optm_Cost_Mtx(i+1=1,:,:)
% jnxt_minCost is the choice of path associated with this
Optm_Path_Mtx(1,1,1) = jnxt_minCost; % choose to store this in first element
Optm_Cost_Mtx(1,1,1) = minCost; 

% clean up and exit (no error)

if run_getOptimalPathMatrix_mex
% note: Optm_Path_Mtx_mex not implemented (due to 1-based vs 0-based
% counting), so just rely on Optm_Cost_Mtx_mex comparison
for j1=1:6
    if(norm(Optm_Cost_Mtx(:,:,j1)-Optm_Cost_Mtx_mex(:,:,j1))>1e-12)
        asdf=5;
        error('getOptimalPathMatrix_mex not working');
    end
end
end

end

function bool=PathIsObstructed(grid_isOccupiedBy,i,j,jnxt,UNOCCUPIED)
% returns true if path from (i,j) to (i+1,jnxt) is obstructed.
% True if any of the points in the box defined by its diagonal corners
% at (i,j) and (i+1,jnxt) is occupied.


if j<=jnxt
    for jj=j:jnxt
        if grid_isOccupiedBy(i,jj) ~= UNOCCUPIED || ...
           grid_isOccupiedBy(i+1,jj) ~= UNOCCUPIED
            bool = true; return;
        end
    end
else
    for jj=j:-1:jnxt
        if grid_isOccupiedBy(i,jj) ~= UNOCCUPIED || ...
           grid_isOccupiedBy(i+1,jj) ~= UNOCCUPIED
            bool = true; return;
        end
    end
end

bool = false;

end

function Cost_AngleChange=getCost_AngleChange(grid_x,grid_y,X,i,j,j1,jnxt,DesignParameters,ModelParams) %#ok<INUSL>
% cost of angle change in going from (i-1,j1) to (i,j) to (i+1,jnxt)
% i goes from 0 to N-1
%
% implemented soft boundaries for angle change at value of <maxanglechange>

% nx=ModelParams.nx;
stateindex_x=ModelParams.stateindex_x;
stateindex_y=ModelParams.stateindex_y;
stateindex_phi=ModelParams.stateindex_phi;
% stateindex_v=ModelParams.stateindex_v;

% anglecost_quadscaling = 100;
% for abs(dtheta)>maxanglechange, have:
% Cost_AngleChange = maxanglechange+anglecost_quadscaling*(abs(dtheta)-maxanglechange)^2;
% with:
%     maxanglechange = 20;
%     anglecost_quadscaling = 100; 
%     DEAD_END_COST = 1000;
% at about 23 degrees, cost is comparable to dead-end cost

if i==0
%     dx1 = X(nx+stateindex_x) - X(stateindex_x);
%     dy1 = X(nx+stateindex_y) - X(stateindex_y);
%     
%     dx2 = grid_x(1,jnxt) - X(stateindex_x);
%     dy2 = grid_y(1,jnxt) - X(stateindex_y);
    phi0=X(stateindex_phi);
    
    phi1=atan2(grid_y(1,jnxt)-X(stateindex_y),grid_x(1,jnxt)-X(stateindex_x));
    
    dphi=phi1-phi0; % radians
else

    if i==1
        dx1 = grid_x(1,j) - X(stateindex_x);
        dy1 = grid_y(1,j) - X(stateindex_y);
    else
        dx1 = grid_x(i,j) - grid_x(i-1,j1);
        dy1 = grid_y(i,j) - grid_y(i-1,j1);
    end

    dx2 = grid_x(i+1,jnxt) - grid_x(i,j);
    dy2 = grid_y(i+1,jnxt) - grid_y(i,j);
    
    dotproduct = dx1*dx2 + dy1*dy2;
    denom = sqrt( (dx1*dx1+dy1*dy1) * (dx2*dx2+dy2*dy2) );
    
    ratio = dotproduct/denom;
    
    if(ratio>1),ratio=1;
    elseif(ratio<-1),ratio=-1;
    end

    dphi = acos( ratio ); % acos ranges from 0 to pi rad
    % dphi = acosd( dotproduct/denom ); % acosd ranges from 0 to 180 deg
    
%     if(i==12 && j==4 && j1==1)
%         fprintf('d1: (%g , %g)\n',dx1,dy1);
%         fprintf('d2: (%g , %g)\n',dx2,dy2);
%         fprintf('ratio = %g / %g = %g\n',dotproduct,denom,ratio);
%         fprintf('dphi = %g\n',dphi);
%     end
end

    
Cost_AngleChange = dphi; % radians
return;



% calculate maximum allowed angle change
deltamax=20*pi/180; %[rad] (defined as input constraint)
speed=X(nx*i+stateindex_v);
try
    maxanglechange = deltamax*speed/ModelParams.L*DesignParameters.Ts;
catch ME
    asdf=5;
    error('maxanglechange formula invalid');
    maxanglechange = 15*pi/180; % [rad]
end

% should under no circumstance choose an angle greater than maxanglechange
if abs(dtheta)>maxanglechange
    %Cost_AngleChange=maxanglechange+anglecost_quadscaling*(abs(dphi)-maxanglechange)^2;
    Cost_AngleChange = DesignParameters.Cost_maxanglechange; 
else
    Cost_AngleChange = dphi;
end

end

function Cost_PathLength=getCost_PathLength(grid_x,grid_y,X,i,j,jnxt,stateindex_x,stateindex_y)
% cost of path length in going from (i,j) to (i+1,jnxt)
% i goes from 0 to N-1

if i==0
    dx = grid_x(i+1,jnxt) - X(stateindex_x);
    dy = grid_y(i+1,jnxt) - X(stateindex_y);
else
    dx = grid_x(i+1,jnxt) - grid_x(i,j);
    dy = grid_y(i+1,jnxt) - grid_y(i,j);
end

d = sqrt(dx^2+dy^2);

Cost_PathLength = d;

end

function Cost_PathDeviation=getCost_PathDeviation(priorPath,i,jnxt)
% cost of deviation from previously planned path when at (i+1,jnxt)
% i goes from 0 to N-1

dj = abs( priorPath(i+1) - jnxt );

% cast to double
dj = cast(dj, 'double');

Cost_PathDeviation = dj*100;


end

function OptmPath=TraceOptimalPath(Optm_Path_Mtx,DesignParameters,Constants)
% Optm_Path_Mtx is [N]by[grid_width]by[grid_width]
% OptmPath is [N]by[1]

N=DesignParameters.N;
DEAD_END_PATH = Constants.DEAD_END_PATH;

OptmPath = zeros(N,1); % initialize.
% OptmPath = cast(OptmPath, 'int32'); % cast to int

OptmPath(1) = Optm_Path_Mtx(1,1,1);
OptmPath(2) = Optm_Path_Mtx(2,OptmPath(1),1);

for i1=3:N
    % i1 represents i+1

    OptmPath(i1) = Optm_Path_Mtx(i1,OptmPath(i1-1),OptmPath(i1-2));
    if OptmPath(i1)==DEAD_END_PATH
        % exit function. Find_Nearest_Borders() will take care of the rest.
        disp('BLAAAAAAAAAAAAAAAAAAAAAAAAAAAAH');
        return;
    end

end

end

function [new_b last_border]=Find_Nearest_Borders(X_all,grid_x,grid_y,new_b_default,grid_isOccupiedBy,OptmPath,DesignParameters,ModelParams,Constants)

% new_b_default is [Lx, Ly, Rx, Ry]
% note: have made assumption that state indices for x,y,psi are 1,2,3

run_Find_Nearest_Borders_mex=0;
if run_Find_Nearest_Borders_mex
[new_b_mex last_border_mex]=Find_Nearest_Borders_mex(X_all,grid_x,grid_y,new_b_default,grid_isOccupiedBy,OptmPath);
end

grid_width=DesignParameters.grid_width;
nx=ModelParams.nx;
N=DesignParameters.N;
carlength=DesignParameters.carlength;
carwidth=DesignParameters.carwidth;
UNOCCUPIED=Constants.UNOCCUPIED;
DEAD_END_PATH = Constants.DEAD_END_PATH;
LAST_BORDER_XMAX=Constants.LAST_BORDER_XMAX;

new_b_L = zeros(N,2);
new_b_R = zeros(N,2);

% for each row
for i=1:N
    j_optm = OptmPath(i);
    
    if j_optm==DEAD_END_PATH
        
        % (i,j_optm) is invalid for any choice of j_optm
        % note that (i-1,j_prev) is still a valid position
        if i==1
            % cannot use previous borders, what to do?
            error('no valid first point');
        end
        
        j_prev=OptmPath(i-1);
        
        % if not first point, use borders from previous point for rest of
        % horizon
        for ii=i:N
            new_b_L(ii,:)=new_b_L(i-1,:);
            new_b_R(ii,:)=new_b_R(i-1,:);
        end
        new_b = [new_b_L new_b_R];
        
        % set the last border
        P1=[grid_x(i-1,j_prev) grid_y(i-1,j_prev)];
        P2=[grid_x(  i,j_prev) grid_y(  i,j_prev)];
        last_border=calculate_last_border(P1,P2);
        
%         figure(5);hold on;
%         plot(grid_x,grid_y,'go');
%         plot([P1(1) P2(1)],[P1(2) P2(2)],'r*');
%         xx=[-2,2];
%         yy=(last_border(3)-last_border(1)*xx)/last_border(2);
%         plot(xx,yy,'b');
%         axis([-2 2 -2 2]);
        
        return; % exit from function
    end
    
    % x,y coordinates of left and right track borders
    Lx=new_b_default(i,1);
    Ly=new_b_default(i,2);
    Rx=new_b_default(i,3);
    Ry=new_b_default(i,4);
    
    % identify which is first point to left that is occupied
    for j=j_optm:-1:1
        k = grid_isOccupiedBy(i,j); % car index
        if k~=UNOCCUPIED
            % this is the first occupied point to left of car
            
            % get states of car
            xypsi = X_all(i*nx+1:i*nx+3,k);

            
            % find closest intersect point to right track border
            new_b_L(i,:) = FindClosestIntersect(Rx,Ry,Lx,Ly,xypsi,carlength,carwidth);

            break;
        elseif j==1
            % has hit track border
            new_b_L(i,:) = new_b_default(i,1:2);
        end
    end
    
    % identify which is first point to right that is occupied
    for j=j_optm:grid_width
        k = grid_isOccupiedBy(i,j); % car index
        if k~=UNOCCUPIED
            % this is the first occupied point to right of car
            
            % get states of car
            xypsi = X_all(i*nx+1:i*nx+3,k);
                        
            % find closest intersect point to left track border
            new_b_R(i,:) = FindClosestIntersect(Lx,Ly,Rx,Ry,xypsi,carlength,carwidth);
                        
            break;
        elseif j==grid_width
            % has hit track border
            new_b_R(i,:) = new_b_default(i,3:4);
        end
    end
end

new_b = [new_b_L new_b_R];

% set the last border as a "trivial" last border
last_border=[1;0;LAST_BORDER_XMAX]; % this constraint means: x<=LAST_BORDER_XMAX

if run_Find_Nearest_Borders_mex
if((norm(new_b - new_b_mex)>1e-12) || norm(last_border-last_border_mex)>1e-12)
   asdf=5;
   error('mex function Find_Nearest_Borders_mex not working.');
end
end

end

function last_border=calculate_last_border(P1,P2)
% Calculates the half-plane constraint: [a1 b2][x;y]<=b
% last_border=[a1;a2;b]
%
% The half-plane is defined by the line passing through the midpoints
% between P1 and P2 and perpendicular to the line through P1 and P2.
%
% choose your ordering of P1 and P2 s.t. P1 will fulfill the constraint, 
% while P2 does not.

x1=P1(1);
x2=P2(1);
y1=P1(2);
y2=P2(2);

a1=x2-x1;
a2=y2-y1;
b=0.5*(-x1^2-y1^2+x2^2+y2^2);

% check if direction of half-plane is correct
if a1*x1+a2*y1<=b
    last_border=[a1;a2;b];
else
    % change the direction of the half-plane
    last_border=[-a1;-a2;-b];
end

end


function intersect = FindClosestIntersect(xA,yA,xB,yB,xypsi,carlength,carwidth)
% finds the intersection between the 4 sides of the car and the line AB.
% in general there are up to two: this finds the intersection closest to A.
%
% note: assumed other car could be backwards, hence need to check all 4 sides
% 
% http://paulbourke.net/geometry/lineline2d/

center=xypsi(1:2);
psi=xypsi(3);

cospsi=cos(psi);
sinpsi=sin(psi);

R = [cospsi -sinpsi;sinpsi cospsi]; % rotation matrix

% first calculate position of 4 corners
% corners is [2]by[4]
l=carlength/2; w=carwidth/2;
corners(:,1) = center + R*[ l; w]; %FL
corners(:,2) = center + R*[ l;-w]; %FR
corners(:,3) = center + R*[-l;-w]; %BR
corners(:,4) = center + R*[-l; w]; %BL

% fAB is the normalized fraction from point A to point B that the intersect
% lies on.
fABinit=2.0;
fABmin=fABinit; % initialize (need only be larger than 1)


% for each adjacent pair of corners, find intersect with line AB
for i=1:4
    if i==4
        j=1;
    else
        j=i+1;
    end
    xi = corners(1,i);
    yi = corners(2,i);
    xj = corners(1,j);
    yj = corners(2,j);
    
    % http://paulbourke.net/geometry/lineline2d/
%     numerij = (xB-xA)*(yA-yi) - (yB-yA)*(xA-xi);
%     denom = (yj-yi)*(xB-xA) - (xj-xi)*(yB-yA);
    
    numerij = (yA-yB)*(xi-xA) + (xB-xA)*(yi-yA);
    denom = (xB-xA)*(yi-yj) - (xi-xj)*(yB-yA);
    
    if denom==0
        %two lines are parallel
        continue;
    end
    
    % normalized distance between point i and point j where intersect lies
    fij = numerij/denom; 
    
    if ~( fij<=1 && fij>=0 )
        % intersect does not lie between point i and j
        continue;
    end
    
    % intersect lies between i,j so see how close it is to A
%     numerAB = (xj-xi)*(yA-yi) - (yj-yi)*(xA-xi);

    numerAB = (yi-yj)*(xi-xA) + (xj-xi)*(yi-yA);
    
    fAB = numerAB/denom;
    
    if ~( fAB<=1 && fAB>=0 )
        % intersect does not lie between point A and B
        continue;
    end
    
    if fAB<fABmin
        fABmin=fAB;
    end
end

if fABmin==fABinit,error('no intersect between line and polygon'),end
intersect = [xA;yA] + fABmin*[xB-xA;yB-yA];
intersect=intersect'; % return a row vector

% plotaid_FindClosestIntersect(xA,yA,xB,yB,corners,intersect)

end

function [in] = inpoly(points,vertices)

% determinse whether a point is inside the polytop given by the vertices
n = length(points);
CorrectSum = sum(convhull(vertices(:,1),vertices(:,2)));
in = false(n,1);
for i = 1:n
    
	tmp_vertices = [vertices;points(i,:)];
    if sum(convhull(tmp_vertices(:,1),tmp_vertices(:,2))) == CorrectSum
        in(i) = 1;
    else
        in(i) = 0;
    end
end

end
