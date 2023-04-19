function generateArcLengthSpline(x,y)
    cleanPath = outlierRemoval(x, y);
    fitSpline(cleanPath.x, cleanPath.y);
end

function resampledPath = outlierRemoval(xOriginal, yOriginal)
    % remove points which are not at all equally spaced, to avoid fitting problems

  % compute mean distance between points and then process the points such that points
  % are not closer than 0.75 the mean distance

  nPoints = size(xOriginal);

  % initialize with zero
  resampledPath.x = zeros(nPoints);
  resampledPath.y = zeros(nPoints);

  % compute distance between points in X-Y data
  distVec = zeros(nPoints - 1);
  for i = 1:(nPoints - 1)
    dx = xOriginal(i + 1) - xOriginal(i);
    dy = yOriginal(i + 1) - yOriginal(i);
    distVec(i) = sqrt(dx * dx + dy * dy);
  end
  % compute mean distance between points
  meanDist = sum(distVec) / (nPoints - 1);

  % compute the new points
  % start point is the original start point
  k = 0;
  resampledPath.x(k) = xOriginal(k);
  resampledPath.y(k) = yOriginal(k);
  k = k + 1;
  for i = 2:(nPoints - 1)
    % compute distance between currently checked point and the one last added to the new X-Y path
    dx = xOriginal(i) - xOriginal(j);
    dy = yOriginal(i) - yOriginal(j);
    dist = sqrt(dx * dx + dy * dy);
    % if this distance is smaller than 0.7 the mean distance add this point to the new X-Y path
    if dist >= 0.7 * meanDist
      resampledPath.x(k) = xOriginal(i);
      resampledPath.y(k) = yOriginal(i);
      k = k+1;
      j = i;
    end
  end
  % always add the last point
  resampledPath.x(k) = xOriginal(nPoints - 1);
  resampledPath.y(k) = yOriginal(nPoints - 1);
  k = k+1;

  % set the new X-Y data
  %    setData(X.head(k),Y.head(k));
  resampledPath.x = reshape(resampledPath.x, [1,k]);
  resampledPath.y = reshape(resampledPath.y, [1,k]);
end

function fitSpline(x,y)
  % successively fit spline -> re-sample path -> compute arc length
  % temporary spline class only used for fitting
  sApproximation = compArcLength(x, y);
  %    std::cout << sApproximation << std::endl;
  totalArcLength = sApproximation(sApproximation.size() - 1);
  
  % 1. spline fit
  firstSplineX.genSpline(sApproximation, x, false);
  firstSplineY.genSpline(sApproximation, y, false);
  % 1. re-sample
  firstRefinedPath = resamplePath(firstSplineX, firstSplineY, totalArcLength);
  sApproximation = compArcLength(firstRefinedPath.x, firstRefinedPath.y);

  totalArcLength = sApproximation(sApproximation.size() - 1);
  %///////////////////////////////////////////
  % 2. spline fit
  secondSplineX.genSpline(sApproximation, firstRefinedPath.x, false);
  secondSplineY.genSpline(sApproximation, firstRefinedPath.y, false);
  % 2. re-sample
  secondRefinedPath = resamplePath(secondSplineX, secondSplineY, totalArcLength);
  %///////////////////////////////////////////
  setRegularData(secondRefinedPath.x, secondRefinedPath.y, secondRefinedPath.s);
  %    setData(secondRefinedPath.X,secondRefinedPath.Y);
  % Final spline fit with fixed Delta_s
  d_splineX.genSpline(d_pathData.s, d_pathData.x, true);
  d_splineY.genSpline(d_pathData.s, d_pathData.y, true);
end

function s = compArcLength(xIn,yIn)
  % compute arc length based on straight line distance between the data points
  nPoints = size(xIn);
  % initailize s as zero
  s = zeros(nPoints);
  %    std::cout << xIn << std::endl;
  for i = 1:(nPoints - 1)
    dx = xIn(i + 1) - xIn(i);
    dy = yIn(i + 1) - yIn(i);
    dist = sqrt(dx * dx + dy * dy);  % dist is straight line distance between points
    s(i + 1) = s(i) + dist;               % s is cumulative sum of dist
  end
  %    std::cout << sIn << std::endl;
end