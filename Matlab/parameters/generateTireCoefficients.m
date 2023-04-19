function tire_coefficients = generateTireCoefficients()
%GENERATETIRECOEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here
    tire_coefficients.scaling_coefficients = generateScalingCoefficients();
    tire_coefficients.longitudinal_coefficients = generateLongitudinalCoefficients();
    tire_coefficients.lateral_coefficients = generateLateralCoefficients();
    tire_coefficients.rolling_coefficients = generateRollingCoefficients();
end

function scaling_coefficients = generateScalingCoefficients()
%GENERATESCALINGCOEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here
  scaling_coefficients.lfzo = 1.0;
  scaling_coefficients.lcx = 1.0;
  scaling_coefficients.lmux = 1.0;
  scaling_coefficients.lex = 1.0;
  scaling_coefficients.lkx = 1.0;
  scaling_coefficients.lhx = 1.0;
  scaling_coefficients.lvx = 1.0;
  scaling_coefficients.lgax = 1.0;
  scaling_coefficients.lcy = 1.0;
  scaling_coefficients.lmuy = 1.0;
  scaling_coefficients.ley = 1.0;
  scaling_coefficients.lky = 1.0;
  scaling_coefficients.lhy = 1.0;
  scaling_coefficients.lvy = 1.0;
  scaling_coefficients.lgay = 1.0;
  scaling_coefficients.ltr = 1.0;
  scaling_coefficients.lres = 0.0;
  scaling_coefficients.lgaz = 1.0;
  scaling_coefficients.lxal = 1.0;
  scaling_coefficients.lyka = 1.0;
  scaling_coefficients.lvyka = 1.0;
  scaling_coefficients.ls = 1.0;
  scaling_coefficients.lsgkp = 1.0;
  scaling_coefficients.lsgal = 1.0;
  scaling_coefficients.lgyr = 1.0;
  scaling_coefficients.lmx = 1.0;
  scaling_coefficients.lvmx = 1.0;
  scaling_coefficients.lmy = 1.0;
end

function longitudinal_coefficients = generateLongitudinalCoefficients()
%GENERATELONGITUDINALCOEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here
  longitudinal_coefficients.pcx1 = +1.786E+000;
  longitudinal_coefficients.pdx1 = +2.850E+000;
  longitudinal_coefficients.pdx2 = -3.540E-001;
  longitudinal_coefficients.pdx3 = +1.223E+001;
  longitudinal_coefficients.pex1 = +8.710E-001;
  longitudinal_coefficients.pex2 = -3.800E-002;
  longitudinal_coefficients.pex3 = +0.000E+000;
  longitudinal_coefficients.pex4 = +7.100E-002;
  longitudinal_coefficients.pkx1 = +8.125E+001;
  longitudinal_coefficients.pkx2 = -2.025E+001;
  longitudinal_coefficients.pkx3 = +5.000E-001;
  longitudinal_coefficients.phx1 = +0.000E+000;
  longitudinal_coefficients.phx2 = +0.000E+000;
  longitudinal_coefficients.pvx1 = +0.000E+000;
  longitudinal_coefficients.pvx2 = +0.000E+000;
  longitudinal_coefficients.rbx1 = +2.372E+001;
  longitudinal_coefficients.rbx2 = +2.597E+001;
  longitudinal_coefficients.rcx1 = +7.495E-001;
  longitudinal_coefficients.rex1 = -4.759E-001;
  longitudinal_coefficients.rex2 = +8.109E-001;
  longitudinal_coefficients.rhx1 = +0.000E+000;
end

function lateral_coefficients = generateLateralCoefficients()
%GENERATELATTERALCOEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here
  lateral_coefficients.pcy1 = +1.434E+000;
  lateral_coefficients.pdy1 = +2.716E+000;
  lateral_coefficients.pdy2 = -5.444E-001;
  lateral_coefficients.pdy3 = +5.190E+000;
  lateral_coefficients.pey1 = -4.869E-001;
  lateral_coefficients.pey2 = -1.487E+000;
  lateral_coefficients.pey3 = +6.282E-002;
  lateral_coefficients.pey4 = +1.154E+000;
  lateral_coefficients.pky1 = -5.322E+001;
  lateral_coefficients.pky2 = +2.060E+000;
  lateral_coefficients.pky3 = +8.336E-001;
  lateral_coefficients.phy1 = +0.000e+000;
  lateral_coefficients.phy2 = +0.000e+000;
  lateral_coefficients.phy3 = -2.030E-002;
  lateral_coefficients.pvy1 = +0.000e+000;
  lateral_coefficients.pvy2 = +0.000e+000;
  lateral_coefficients.pvy3 = -2.713E+000;
  lateral_coefficients.pvy4 = -1.517E+000;
  lateral_coefficients.rby1 = +2.033e+001;
  lateral_coefficients.rby2 = +8.152e+000;
  lateral_coefficients.rby3 = -1.243e-002;
  lateral_coefficients.rcy1 = +9.317e-001;
  lateral_coefficients.rey1 = -3.982e-004;
  lateral_coefficients.rey2 = +3.077e-001;
  lateral_coefficients.rhy1 = +0.000e+000;
  lateral_coefficients.rhy2 = +0.000e+000;
  lateral_coefficients.rvy1 = +0.000e+000;
  lateral_coefficients.rvy2 = +0.000e+000;
  lateral_coefficients.rvy3 = +0.000e+000;
  lateral_coefficients.rvy4 = +0.000e+000;
  lateral_coefficients.rvy5 = +0.000e+000;
  lateral_coefficients.rvy6 = +0.000e+000;
end

function rolling_coefficients = generateRollingCoefficients()
%GENERATEROLLINGCOEFFICIENTS Summary of this function goes here
%   Detailed explanation goes here
  rolling_coefficients.qsy1 = -2.367E-002;
  rolling_coefficients.qsy2 = +0.000e+000;
  rolling_coefficients.qsy3 = +0.000e+000;
  rolling_coefficients.qsy4 = +0.000e+000;
end

