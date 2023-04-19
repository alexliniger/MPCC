classdef SplineData
    properties (Access = public)
        xData
        yData
        nPoints
        isRegular
        deltaX
        xMap
    end

    methods
        function obj = SplineData()
            obj.xMap = containers.Map('KeyType','double','ValueType','int64');
        end
    end
end
