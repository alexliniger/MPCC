classdef Bounds
    %BOUNDS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        d_uBoundsX
        d_lBoundsX
        d_uBoundsU
        d_lBoundsU
        d_uBoundsS
        d_lBoundsS
    end
    
    methods (Access = public)
        function obj = Bounds(bounds)
            obj.d_uBoundsX = zeros(11,1);
            obj.d_lBoundsX = zeros(11,1);
            obj.d_uBoundsU = zeros(4,1);
            obj.d_lBoundsU = zeros(4,1);
            obj.d_uBoundsS = zeros(2,1);
            obj.d_lBoundsS = zeros(2,1);

            obj.d_lBoundsX(1) = bounds.lowerStateBounds.xL;
            obj.d_lBoundsX(2) = bounds.lowerStateBounds.yL;
            obj.d_lBoundsX(3) = bounds.lowerStateBounds.phiL;
            obj.d_lBoundsX(4) = bounds.lowerStateBounds.vxL;
            obj.d_lBoundsX(5) = bounds.lowerStateBounds.vyL;
            obj.d_lBoundsX(6) = bounds.lowerStateBounds.rL;
            obj.d_lBoundsX(7) = bounds.lowerStateBounds.sL;
            obj.d_lBoundsX(8) = bounds.lowerStateBounds.throttleL;
            obj.d_lBoundsX(9) = bounds.lowerStateBounds.steeringAngleL;
            obj.d_lBoundsX(10) = bounds.lowerStateBounds.brakesL;
            obj.d_lBoundsX(11) = bounds.lowerStateBounds.vsL;

            obj.d_uBoundsX(1) = bounds.upperStateBounds.xU;
            obj.d_uBoundsX(2) = bounds.upperStateBounds.yU;
            obj.d_uBoundsX(3) = bounds.upperStateBounds.phiU;
            obj.d_uBoundsX(4) = bounds.upperStateBounds.vxU;
            obj.d_uBoundsX(5) = bounds.upperStateBounds.vyU;
            obj.d_uBoundsX(6) = bounds.upperStateBounds.rU;
            obj.d_uBoundsX(7) = bounds.upperStateBounds.sU;
            obj.d_uBoundsX(8) = bounds.upperStateBounds.throttleU;
            obj.d_uBoundsX(9) = bounds.upperStateBounds.steeringAngleU;
            obj.d_uBoundsX(10) = bounds.upperStateBounds.brakesU;
            obj.d_uBoundsX(11) = bounds.upperStateBounds.vsU;

            obj.d_lBoundsU(1) = bounds.lowerInputBounds.dThrottleL;
            obj.d_lBoundsU(2) = bounds.lowerInputBounds.dSteeringAngleL;
            obj.d_lBoundsU(3) = bounds.lowerInputBounds.dBrakesL;
            obj.d_lBoundsU(4) = bounds.lowerInputBounds.dVsL;

            obj.d_uBoundsU(1) = bounds.upperInputBounds.dThrottleU;
            obj.d_uBoundsU(2) = bounds.upperInputBounds.dSteeringAngleU;
            obj.d_uBoundsU(3) = bounds.upperInputBounds.dBrakesU;
            obj.d_uBoundsU(4) = bounds.upperInputBounds.dVsU;
        end

        function boundsLX = getBoundsLX(obj,x)
            boundsLX = obj.d_lBoundsX - x;
        end
        
        function boundsUX = getBoundsUX(obj,x) 
            boundsUX = obj.d_uBoundsX - x;
        end
        
        function boundsLU = getBoundsLU(obj,u)
            boundsLU = obj.d_lBoundsU - u;
        end
        
        function boundsUU = getBoundsUU(obj,u) 
            boundsUU = obj.d_uBoundsU - u;
        end
        
        function boundsLS = getBoundsLS(obj)
            boundsLS = obj.d_lBoundsS;
        end
        
        function boundsUS = getBoundsUS(obj)
            boundsUS = obj.d_uBoundsS;
        end
    end
end

