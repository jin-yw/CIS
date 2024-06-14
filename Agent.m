classdef Agent < handle
    % The agent model
    
    properties
        Range;
        MeasurementModel;
        Pd;
    end
    
    methods
        function td = Agent(range, measurementmodel, pd)
            td.Range = range;
            td.MeasurementModel = measurementmodel;
            td.Pd = pd;
        end
        function [pd] = getPd(this, x)
            z = this.MeasurementModel.estimateNewState(x);
            zLen = length(z);
            flag = true;
            for idx = 1:zLen
                if z(idx) < this.Range(idx,1) || z(idx) > this.Range(idx,2)
                    flag = false;
                    break;
                end
            end
            if flag
                pd = this.Pd;
            else
                pd = 0;
            end
        end
    end
    
end

