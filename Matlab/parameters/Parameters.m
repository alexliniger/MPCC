classdef Parameters
    %PARAMETERS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        bounds
        mpcModel
        costs
        car
        tire
        normalization
        config
    end
    
    methods
        function obj = Parameters()
            %load bounds
            fname = 'bounds.json';
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            obj.bounds = jsondecode(str);

            %load mpcmodelparameters
            fname = 'model.json';
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            obj.mpcModel = jsondecode(str);

            %load costs
            fname = 'cost.json';
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            obj.costs = jsondecode(str);

            %load car parameters
            fname = 'car.json';
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            obj.car = jsondecode(str);

            %load tire coefficients
            fname = 'tire.json';
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            obj.tire = jsondecode(str);

            %load noermalization parameters
            fname = 'normalization.json';
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            obj.normalization = jsondecode(str);

            %load config
            fname = 'config.json';
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            obj.config = jsondecode(str);
        end
    end
end

