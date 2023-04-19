function spline = genSpline(xIn, yIn, isRegular)
{
  % given x and y data generate spline
  % special case for regular or irregular spaced data points in x
  % if regular the spacing in x is given by deltaX

  % store data in data struct
  if isRegular
    deltaX = xIn(1) - xIn(0);
    setRegularData(xIn, yIn, deltaX);
  else 
    setData(xIn, yIn);
  end
  % given data compute spline parameters

  % TODO if success is false call exception

  % bool success = compSplineParams();

  spline = compSplineParams();
}

