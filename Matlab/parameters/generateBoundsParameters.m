function bounds_parameters = generateBoundsParameters()
    %GENERATEBOUNDSPARAMETERS Summary of this function goes here
    %   Detailed explanation goes here
    bounds_parameters.l_state_bounds = generateLowerStateBounds();
    bounds_parameters.u_state_bounds = generateUpperStateBounds();
    
    bounds_parameters.l_input_bounds = generateLowerInputBounds();
    bounds_parameters.u_input_bounds = generateUpperInputBounds();
end

function l_state_bounds = generateLowerStateBounds()

    l_state_bounds.xL = -3000.0;
    l_state_bounds.yL = -3000.0;
    l_state_bounds.phiL = -10.0;
    l_state_bounds.vxL = 0.0;
    l_state_bounds.vyL = -10.0;
    l_state_bounds.rL = -1.0;
    l_state_bounds.sL = -500.0;
    l_state_bounds.throttleL = 0.0;
    l_state_bounds.steeringAngleL = -1.0;
    l_state_bounds.brakesL = 0.0;
    l_state_bounds.vsL = 0.0;

end

function u_state_bounds = generateUpperStateBounds()

    u_state_bounds.xU = 3000.0;
    u_state_bounds.yU = 3000.0;
    u_state_bounds.phiU = 10.0;
    u_state_bounds.vxU = 30.0;
    u_state_bounds.vyU = 10.0;
    u_state_bounds.rU = 2.0;
    u_state_bounds.sU = 2500.0;
    u_state_bounds.throttleU = 1.0;
    u_state_bounds.steeringAngleU = 1.0;
    u_state_bounds.brakesU = 1.0;
    u_state_bounds.vsU = 50.0;

end

function l_input_bounds = generateLowerInputBounds()

    l_input_bounds.dThrottleL = -1.0;
    l_input_bounds.dSteeringAngleL = -1.0;
    l_input_bounds.dBrakesL = -1.0;
    l_input_bounds.dVsL = -50.0;

end

function u_input_bounds = generateUpperInputBounds()

    u_input_bounds.dThrottleU = 1.0;
    u_input_bounds.dSteeringAngleU = 1.0;
    u_input_bounds.dBrakesU = 1.0;
    u_input_bounds.dVsU = 50.0;

end





