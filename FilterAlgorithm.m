classdef FilterAlgorithm < handle
    % 算法类的基类

    properties
        KinematicModel;         % 运动模型
        MeasurementModel;       % 测量模型
    end
    
    methods
        % 设置时间的函数，步长一般用于运动模型中
        function setT(obj,t)
            if t > 0
                obj.KinematicModel.setT(t);
%                 obj.MeasurementModel.setT(t);
            end
        end
    end
    
end

