function normalization_parameters = generateNormalizationParameters()
%GENERATENORMALIZATIONPARAMETERS Summary of this function goes here
%   Detailed explanation goes here
    x = 1.0;
    y = 1.0;
    yaw = 1.0;
    vx = 1.0;
    vy = 1.0;
    r = 1.0;
    s = 1.0;
    throttle = 1.0;
    steeringAngle = 1.0;
    brakes = 1.0;
    vs = 1.0;
    dThrottle = 1.0;
    dSteeringAngle = 1.0;
    dBrakes = 1.0;
    dVs = 1.0;

    norm_state_vec = [x, y, yaw, vx, vy, r, s, throttle, steeringAngle, brakes, vs];
    inv_norm_state_vec = [1.0 / x, 1.0 / y, 1.0 / yaw, 1.0 / vx, 1.0 / vy, 1.0 / r, 1.0 / s, 1.0 / throttle, 1.0 / steeringAngle, 1.0 / brakes, 1.0 / vs];
    norm_input_vec = [dThrottle, dSteeringAngle, dBrakes, dVs];
    inv_norm_input_vec = [1.0 / dThrottle, 1.0 / dSteeringAngle, 1.0 / dBrakes, 1.0 / dVs];

    normalization_parameters.tX = diag(norm_state_vec);
    normalization_parameters.tXInv = diag(inv_norm_state_vec);

    normalization_parameters.tU = diag(norm_input_vec);
    normalization_parameters.tUInv = diag(inv_norm_input_vec);

    normalization_parameters.tS = eye(2);
    normalization_parameters.tSInv = eye(2);
end

