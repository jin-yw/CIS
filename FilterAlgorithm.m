classdef FilterAlgorithm < handle
    % �㷨��Ļ���

    properties
        KinematicModel;         % �˶�ģ��
        MeasurementModel;       % ����ģ��
    end
    
    methods
        % ����ʱ��ĺ���������һ�������˶�ģ����
        function setT(obj,t)
            if t > 0
                obj.KinematicModel.setT(t);
%                 obj.MeasurementModel.setT(t);
            end
        end
    end
    
end

