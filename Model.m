classdef Model < handle
    % 模型类的基类
    
    properties (SetAccess = protected)
        Std;            % 各个维度上的标准差
        T;              % 仿真过程中，两帧之间的时间间隔
        CovMatrix;      % 协方差矩阵
    end
    
    methods
    end
    
end

