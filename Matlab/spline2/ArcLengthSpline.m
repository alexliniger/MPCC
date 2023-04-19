classdef ArcLengthSpline
    %ARC_LENGTH_SPLINE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        d_pathData
        d_splineX
        d_splineY
        d_parameters
    end
    
    methods (Access = public)
      function obj = ArcLengthSpline(modelParameters)
          obj.d_pathData = PathData();
          obj.d_splineX = CubicSpline();
          obj.d_splineY = CubicSpline();
          obj.d_parameters = modelParameters;
      end

          % X and Y spline used for final spline fit
      function obj = gen2DSpline(obj,x,y)
          % generate 2-D arc length parametrized spline given X-Y data
          % remove outliers, depending on how iregular the points are this can help
          cleanPath = obj.outlierRemoval(x, y);
          % successively fit spline and re-sample
          obj = obj.fitSpline(cleanPath.x, cleanPath.y);
      end

      function sPath = getPostion(obj,s)
          sPath(1) = obj.d_splineX.getPoint(s);
          sPath(2) = obj.d_splineY.getPoint(s);
      end

      function dsPath = getDerivative(obj,s)
          dsPath(1) = obj.d_splineX.getDerivative(s);
          dsPath(2) = obj.d_splineY.getDerivative(s);
      end

      function ddsPath = getSecondDerivative(obj,s)
          ddsPath(1) = obj.d_splineX.getSecondDerivative(s);
          ddsPath(2) = obj.d_splineY.getSecondDerivative(s);
      end

      function length = getLength(obj)
          length = obj.d_pathData.s(obj.d_pathData.nPoints);
      end

      function sGuess = projectOnSpline(obj,x)
          pos(1) = x.x;
          pos(2) = x.y;
          sGuess = x.s;
          posPath = obj.getPostion(sGuess);
        
          sOpt = sGuess;
          dist = norm(pos - posPath);
        
          if dist >= obj.d_parameters.maxDistProj
            % std::cout << "dist too large" << std::endl;
            diffXAll = obj.d_pathData.x -pos(1);
            diffYAll = obj.d_pathData.y - pos(2);
            distSquare = diffXAll.^2 + diffYAll.^2;
            [~, minIndex] = min(distSquare);
            sOpt = obj.d_pathData.s(minIndex);
          end
          sOld = sOpt;
          for i = 1:20
            posPath = obj.getPostion(sOpt);
            dsPath = obj.getDerivative(sOpt);
            ddsPath = obj.getSecondDerivative(sOpt);
            diff = posPath - pos;
            jac = 2.0 * diff(1) * dsPath(1) + 2.0 * diff(2) * dsPath(2);
            hessian = 2.0 * dsPath(1) * dsPath(1) + 2.0 * diff(1) * ddsPath(1)...
                             + 2.0 * dsPath(2) * dsPath(2) + 2.0 * diff(2) * ddsPath(2);
            % Newton method
            sOpt = sOpt - jac / hessian;
            sOpt = obj.unwrapInput(sOpt);

            %        std::cout << std::abs(sOld - sOpt) << std::endl;
            if (abs(sOld - sOpt) <= 1e-5)
                sGuess = sOpt;
                return;
            end
            sOld = sOpt;
          end
          % something is strange if it did not converge within 20 iterations, give back the initial guess
      end

      function path = getPath(obj)
          path = obj.d_pathData;
      end
        
    end

    methods (Access = private)
        function obj = setData(obj,xIn,yIn)
          % set input data if x and y have same length
          % compute arc length based on an piecewise linear approximation
          if size(xIn,2) == size(yIn,2)
            obj.d_pathData.x = xIn;
            obj.d_pathData.y = yIn;
            obj.d_pathData.nPoints = size(xIn,2);
            obj.d_pathData.s = compArcLength(xIn, yIn);
          else
            disp("input data does not have the same length");
          end
        end

        function obj = setRegularData(obj,xIn,yIn,sIn)
          % set final x-y data if x and y have same length
          % x-y points are space such that they are very close to arc length parametrized
          if size(xIn,2) == size(yIn,2)
            obj.d_pathData.x = xIn;
            obj.d_pathData.y = yIn;
            obj.d_pathData.nPoints = size(xIn,2);
            obj.d_pathData.s = sIn;
          else
            disp("input data does not have the same length");
          end
        end


        function s = compArcLength(obj,xIn,yIn)
        % compute arc length based on straight line distance between the data points
          nPoints = size(xIn,2);
          % initailize s as zero
     
          s = zeros(1,nPoints);
          %    std::cout << xIn << std::endl;
          for i = 1:(nPoints - 1)
            dx = xIn(1,i + 1) - xIn(1,i);
            dy = yIn(1,i + 1) - yIn(1,i);
            dist = sqrt(dx * dx + dy * dy);  % dist is straight line distance between points
            s(1,i + 1) = s(1,i) + dist;               % s is cumulative sum of dist
          end
          %    std::cout << sIn << std::endl;
        end

        function resampledPath = resamplePath(obj,initialSplineX,initialSplineY,totalArcLength)
            % re-sample arc length parametrized X-Y spline path with N_spline data points
            % using equidistant arc length values
            % successively re-sample, computing the arc length and then fit the path should
            % result in close to equidistant points w.r.t. arc length
            
            % s -> "arc length" where points should be extracted
            % equilly spaced between 0 and current length of path
            NSpline = 5000;
            resampledPath.nPoints = NSpline;
            resampledPath.s = linspace(0, totalArcLength, NSpline);
            
            % initialize new points as zero
            resampledPath.x = zeros(1,NSpline);
            resampledPath.y = zeros(1,NSpline);
            
            % extract X-Y points
            for i = 1:NSpline
                resampledPath.x(1,i) = initialSplineX.getPoint(resampledPath.s(1,i));
                resampledPath.y(1,i) = initialSplineY.getPoint(resampledPath.s(1,i));
            end
        end

        function resampledPath = outlierRemoval(obj,xOriginal,yOriginal)
            % remove points which are not at all equally spaced, to avoid fitting problems

              % compute mean distance between points and then process the points such that points
              % are not closer than 0.75 the mean distance
              k = 1;  % indecies
              j = 1;
            
              if size(xOriginal,2) ~= size(yOriginal,2)
                disp('error')
              end
              %    std::cout << xOriginal << std::endl;
            
              nPoints = size(xOriginal,2);
            
              % initialize with zero
              resampledPath.x = zeros(1,nPoints);
              resampledPath.y = zeros(1,nPoints);
            
              % compute distance between points in X-Y data
              distVec = zeros(1,nPoints - 1);
              for i = 1:(nPoints - 1)
                dx = xOriginal(1,i + 1) - xOriginal(1,i);
                dy = yOriginal(1,i + 1) - yOriginal(1,i);
                distVec(1,i) = sqrt(dx * dx + dy * dy);
              end
              % compute mean distance between points
              meanDist = sum(distVec) / (nPoints - 1);
            
              % compute the new points
              % start point is the original start point
              resampledPath.x(1,k) = xOriginal(1,k);
              resampledPath.y(1,k) = yOriginal(1,k);
              k = k+1;
              for i = 2:(nPoints - 1)
                % compute distance between currently checked point and the one last added to the new X-Y path
                dx = xOriginal(1,i) - xOriginal(1,j);
                dy = yOriginal(1,i) - yOriginal(1,j);
                dist = sqrt(dx * dx + dy * dy);
                % if this distance is smaller than 0.7 the mean distance add this point to the new X-Y path
                if dist >= 0.7 * meanDist
                  resampledPath.x(1,k) = xOriginal(1,i);
                  resampledPath.y(1,k) = yOriginal(1,i);
                  k = k+1;
                  j = i;
                end
              end
              % always add the last point
              resampledPath.x(1,k) = xOriginal(1,nPoints);
              resampledPath.y(1,k) = yOriginal(1,nPoints);
              k = k+1;
            
              %    std::cout << "not resiszed " << X_new.transpose() << std::endl;
              % set the new X-Y data
              %    setData(X.head(k),Y.head(k));
              resampledPath.x(:,k:1:end) = [];
              resampledPath.y(:,k:1:end) = [];
            
              %    std::cout << "resiszed " << X_new.transpose()  << std::endl;
        end

        function obj = fitSpline(obj,x,y)
              % successively fit spline -> re-sample path -> compute arc length
              % temporary spline class only used for fitting
              sApproximation = obj.compArcLength(x, y);
              %    std::cout << sApproximation << std::endl;
              totalArcLength = sApproximation(size(sApproximation,2));
            
              firstSplineX = CubicSpline();
              firstSplineY = CubicSpline();
              secondSplineX = CubicSpline();
              secondSplineY = CubicSpline();
              % 1. spline fit
              firstSplineX = firstSplineX.genSpline(sApproximation, x, false);
              firstSplineY = firstSplineY.genSpline(sApproximation, y, false);
              % 1. re-sample
              firstRefinedPath = obj.resamplePath(firstSplineX, firstSplineY, totalArcLength);
              sApproximation = obj.compArcLength(firstRefinedPath.x, firstRefinedPath.y);
            
              totalArcLength = sApproximation(size(sApproximation,2));
              %//////////////////////////////////////////
              % 2. spline fit
              secondSplineX = secondSplineX.genSpline(sApproximation, firstRefinedPath.x, false);
              secondSplineY = secondSplineY.genSpline(sApproximation, firstRefinedPath.y, false);
              % 2. re-sample
              secondRefinedPath = obj.resamplePath(secondSplineX, secondSplineY, totalArcLength);
              %//////////////////////////////////////////
              obj = obj.setRegularData(secondRefinedPath.x, secondRefinedPath.y, secondRefinedPath.s);
              %    setData(secondRefinedPath.X,secondRefinedPath.Y);
              % Final spline fit with fixed Delta_s
              obj.d_splineX = obj.d_splineX.genSpline(obj.d_pathData.s, obj.d_pathData.x, true);
              obj.d_splineY = obj.d_splineY.genSpline(obj.d_pathData.s, obj.d_pathData.y, true);
        end

        function input = unwrapInput(obj,x)
            xMax = obj.getLength();
            input = x - xMax * floor(x / xMax);
        end
    end
end

