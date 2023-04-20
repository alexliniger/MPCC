addpath('parameters');

parameters = Parameters();

bounds = Bounds(parameters.bounds);

x = zeros(11,1);
u = zeros(4,1);

disp('LX Bounds:');
disp(bounds.getBoundsLX(x));
disp('UX Bounds:');
disp(bounds.getBoundsUX(x));
disp('LU Bounds:');
disp(bounds.getBoundsLU(u));
disp('UU Bounds:');
disp(bounds.getBoundsUU(u));
disp('LS Bounds:');
disp(bounds.getBoundsLS());
disp('US Bounds:');
disp(bounds.getBoundsUS());