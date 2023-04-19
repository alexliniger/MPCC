function mpc_model_parameters = generateMpcModelParameters()
%GENERATEMPCMODELPARAMETERS Summary of this function goes here
%   Detailed explanation goes here
    mpc_model_parameters.carL = 2.0;
    mpc_model_parameters.carW = 1.0;
    
    mpc_model_parameters.rIn = 2.0;
    mpc_model_parameters.rOut = 2.0;
    
    mpc_model_parameters.maxDistProj = 10.0;
    
    mpc_model_parameters.maxAlpha = 0.15;
    
    mpc_model_parameters.initialVelocity = 0.0;
    mpc_model_parameters.sTrustRegion = 30.0;
    
    mpc_model_parameters.vxZero = 3.0;
end

