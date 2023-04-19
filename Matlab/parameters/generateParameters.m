function [tire,car,bounds_parameters,cost_parameters,normalization_parameters,mpc_model_parameters] = generateParameters()
%GENERATEPARAMETERS Summary of this function goes here
%   Detailed explanation goes here
    tire = generateTireCoefficients();
    car = generateCarModelParameters();
    bounds_parameters = generateBoundsParameters();
    cost_parameters = generateCostParameters();
    normalization_parameters = generateNormalizationParameters();
    mpc_model_parameters = generateMpcModelParameters();
end