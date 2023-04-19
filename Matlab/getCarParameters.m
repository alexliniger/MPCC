function car_parameters = getCarParameters()
%GETCARPARAMETERS Summary of this function goes here
%   Detailed explanation goes here
    car_parameters.m = 300.0;
    car_parameters.g = 9.81;
    car_parameters.iz = 110.0;

    car_parameters.lf = 0.765;
    car_parameters.lr = 0.765;

    car_parameters.fz_nominal = 800.0;
    
    car_parameters.cl = -1.3;
    car_parameters.cd = -1.0;
    car_parameters.cm = 3.0;

    car_parameters.cm1 = 240;
    car_parameters.cbf = 1250.0;
    car_parameters.cbr = 612.0;
    car_parameters.iw = 0.24;
    car_parameters.r_dyn = 0.231;
    car_parameters.gear_ratio = 3.5;
    car_parameters.p_max = 60000.0;
end

