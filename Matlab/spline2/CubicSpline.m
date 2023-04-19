classdef CubicSpline
    %CUBIC_SPLINE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        d_dataSet
        d_splineParameters
        d_splineData
    end

    methods (Access = public)
        function obj = CubicSpline()
            obj.d_splineParameters = SplineParams();
            obj.d_splineData = SplineData();
        end

        function obj = genSpline(obj,xIn,yIn,isRegular)
        % given x and y data generate spline
        % special case for regular or irregular spaced data points in x
        % if regular the spacing in x is given by deltaX
        
        % store data in data struct
          if isRegular
            deltaX = xIn(1,2) - xIn(1,1);
            obj = obj.setRegularData(xIn, yIn, deltaX);
          else
            obj = obj.setData(xIn, yIn);
          end
          % given data compute spline parameters
        
          % TODO if success is false call exception
        
          % bool success = compSplineParams();
        
          obj = obj.compSplineParams();
        end
        
        function point = getPoint(obj, x)
        % evaluation of spline a x
          % wrape input to data -> x data needs start at 0 and contain end point!!!
          x = obj.unwrapInput(x);
          % compute index
          index = obj.getIndex(x);
          % access previous points
          xI = obj.d_splineData.xData(1,index);
          % compute diff to point and it's powers
          dx = x - xI;
          dx2 = dx * dx;
          dx3 = dx * dx2;
          % return spline value y = a + b dx + c dx^2 + d dx^3
          point = obj.d_splineParameters.a(1,index) + obj.d_splineParameters.b(1,index) * dx...
                 + obj.d_splineParameters.c(1,index) * dx2 + obj.d_splineParameters.d(1,index) * dx3;
        end

        function derivative = getDerivative(obj,x)
        % evaluate first derivative of spline
          % identical to noram spline with y' = b + 2 c dx + 3 d dx^2
          x = obj.unwrapInput(x);
          index = obj.getIndex(x);
          xI = obj.d_splineData.xData(index);
        
          dx = x - xI;
          dx2 = dx * dx;
          % y' = b + 2 c dx + 3 d dx^2
          derivative = obj.d_splineParameters.b(index) + 2.0 * obj.d_splineParameters.c(index) * dx...
                 + 3.0 * obj.d_splineParameters.d(index) * dx2;
        end

        function secondDerivative = getSecondDerivative(obj,x)
          % evaluate second derivative of spline
          % identical to noram spline with y' = 2 c + 6 d dx
          x = obj.unwrapInput(x);
          index = obj.getIndex(x);
          xI = obj.d_splineData.xData(index);
        
          dx = x - xI;
          % y' = 2 c + 6 d dx
          secondDerivative = 2.0 * obj.d_splineParameters.c(index) + 6.0 * obj.d_splineParameters.d(index) * dx;
        end
    end

    methods (Access = private)
        function obj = setRegularData(obj,xIn,yIn,deltaX)
            % if x and y have same length, stare given data in spline data struct
              if size(xIn,2) == size(yIn,2)
                obj.d_splineData.xData = xIn;
                obj.d_splineData.yData = yIn;
                obj.d_splineData.nPoints = size(xIn,2);
                obj.d_splineData.isRegular = true;
                obj.d_splineData.deltaX = deltaX;
            
                obj.d_dataSet = true;
              else
                disp('input data does not have the same length');
              end
        end
        
        function obj = setData(obj, xIn,yIn)
            % if x and y have same length, stare given data in spline data struct
              if size(xIn,2) == size(yIn,2)
                obj.d_splineData.xData = xIn;
                obj.d_splineData.yData = yIn;
                obj.d_splineData.nPoints = size(xIn,2);
                obj.d_splineData.isRegular = false;
                obj.d_splineData.deltaX = 0;
                for i = 1:size(xIn,2)
                  obj.d_splineData.xMap(xIn(1,i)) = i;
                end
            
                obj.d_dataSet = true;
              else
                disp('input data does not have the same length');
              end
        end

        function obj = compSplineParams(obj)
            % compute spline parameters parameters
            % code is a replica of the wiki code
              % spline parameters from parameter struct initialized to zero
              obj.d_splineParameters.a = zeros(1,obj.d_splineData.nPoints);
              obj.d_splineParameters.b = zeros(1,obj.d_splineData.nPoints - 1);
              obj.d_splineParameters.c = zeros(1,obj.d_splineData.nPoints);
              obj.d_splineParameters.d = zeros(1,obj.d_splineData.nPoints - 1);
            
              % additional variables used to compute a,b,c,d
              mu = zeros(1,obj.d_splineData.nPoints - 1);
              h = zeros(1,obj.d_splineData.nPoints - 1);
              alpha = zeros(1,obj.d_splineData.nPoints - 1);
              l = zeros(1,obj.d_splineData.nPoints);
              z = zeros(1,obj.d_splineData.nPoints);
            
              % a is equal to y data
              obj.d_splineParameters.a = obj.d_splineData.yData;
              % compute h as diff of x data
              for i = 1:(obj.d_splineData.nPoints - 1)
                h(i) = obj.d_splineData.xData(i + 1) - obj.d_splineData.xData(i);
              end
              % compute alpha
              for i = 2:(obj.d_splineData.nPoints - 1)
                alpha(i) = 3.0 / h(i) * (obj.d_splineParameters.a(i + 1) - obj.d_splineParameters.a(i))...
                           -3.0 / h(i - 1) * (obj.d_splineParameters.a(i) - obj.d_splineParameters.a(i - 1));
              end
            
              % compute l, mu, and z
              l(1) = 1.0;
              mu(1) = 0.0;
              z(1) = 0.0;
              for i = 2:(obj.d_splineData.nPoints - 1)
                l(i) = 2.0 * (obj.d_splineData.xData(i + 1) - obj.d_splineData.xData(i - 1)) - h(i - 1) * mu(i - 1);
                mu(i) = h(i) / l(i);
                z(i) = (alpha(i) - h(i - 1) * z(i - 1)) / l(i);
              end
              %l(obj.d_splineData.nPoints - 1) = 1.0;
              z(obj.d_splineData.nPoints) = 0.0;
            
              % compute b,c,d data given the previous work
              obj.d_splineParameters.c(obj.d_splineData.nPoints) = 0.0;
            
              for i = (obj.d_splineData.nPoints - 1):-1:1
                obj.d_splineParameters.c(i) = z(i) - mu(i) * obj.d_splineParameters.c(i + 1);
                obj.d_splineParameters.b(i) =...
                  (obj.d_splineParameters.a(i + 1) - obj.d_splineParameters.a(i)) / h(i)...
                  -(h(i) * (obj.d_splineParameters.c(i + 1) + 2.0 * obj.d_splineParameters.c(i))) / 3.0;
                obj.d_splineParameters.d(i) =...
                  (obj.d_splineParameters.c(i + 1) - obj.d_splineParameters.c(i)) / (3.0 * h(i));
              end
        end
        
        function index = getIndex(obj,x)
              % given a x value find the closest point in the spline to evaluate it
              % special case if x is regularly space
              % assumes wrapped data!
            
              % if special case of end points
              if x == obj.d_splineData.xData(1,obj.d_splineData.nPoints)
                index = obj.d_splineData.nPoints;
                return;
              end
              % if regular index can be found by rounding
              if obj.d_splineData.isRegular
                index = cast(floor(x / obj.d_splineData.deltaX), "int64");
                return;
              % if irregular index need to be searched
              else
                keyVec = keys(obj.d_splineData.xMap);
                arr = cell2mat(keyVec);
                val = arr(arr > x);
                if isempty(val)
                  index  = -1;
                  return
                else
                  index = obj.d_splineData.xMap(val(1,1)) - 1;
                  return;
                end
             end
        end

        function input = unwrapInput(obj,x)
            xMax = obj.d_splineData.xData(obj.d_splineData.nPoints);
            input = x - xMax * floor(x / xMax);
        end
    end
end

