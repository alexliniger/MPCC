function cost_parameters = generateCostParameters()
%GENERATECOSTPARAMETERS Summary of this function goes here
%   Detailed explanation goes here
    cost_parameters.qC = 0.01;
    cost_parameters.qL = 1000.0;
    cost_parameters.qVs = 0.2;
    
    cost_parameters.qMu = 0.1;
    
    cost_parameters.qR = 1E-4;
    
    cost_parameters.qBeta = 5E0;
    cost_parameters.betaKinCost = 1;
    
    cost_parameters.rThrottle = 1E-3;
    cost_parameters.rSteeringAngle = 1E-9;
    cost_parameters.rBrakes = 1E-3;
    cost_parameters.rVs = 1E-5;
    
    cost_parameters.rdThrottle = 1E-1;
    cost_parameters.rdSteeringAngle = 5E2;
    cost_parameters.rdBrakes = 1E-1;
    cost_parameters.rdVs = 1E-4;
    
    cost_parameters.qCNMult = 1000.0;
    cost_parameters.qRNMult = 10.0;
    
    cost_parameters.scQuadTrack = 100.0;
    cost_parameters.scQuadTire = 100.0;
    cost_parameters.scQuadAlpha = 100.0;
    
    cost_parameters.scLinTrack = 10.0;
    cost_parameters.scLinTire = 10.0;
    cost_parameters.scLinAlpha = 10.0;

end

