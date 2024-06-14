classdef DistributedFusionAlgorithm < handle
    %UNTITLED3 此处显示有关此类的摘要
    %   此处显示详细说明
    
    properties (SetAccess = protected)
        Time;
        THETA;
        L;
    end
    
    methods
        function td = DistributedFusionAlgorithm()
            td.Time = 1;
            td.THETA = 1e-5;
            td.L = 10;
        end
    end
end

