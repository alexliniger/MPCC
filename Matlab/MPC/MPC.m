classdef MPC
    %MPC Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Access = private)
        d_ts;
        d_validInitialGuess;
    
        d_stages;
        d_initialGuess;
        d_optimalSolution;
    
        d_nSqp;
        d_sqpMixing;
        d_nNonSolves;
        d_nNoSolvesSqp;
        d_nReset;
    
        d_modelPtr;
        d_costPtr;
        d_constraintsPtr;
        d_trackPtr;
        d_boundsPtr;
        d_normalizationParametersPtr;
        d_modelParametersPtr;
        d_carParametersPtr;
        d_solverInterface;
    end
    
    methods
        function obj = MPC(inputArg1,inputArg2)
            %MPC Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

