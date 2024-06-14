classdef CVmodel < Model
    % ���ٶ�ģ��
    % ״̬��Ϊ��
    % state = [x, vx, y, vy, z, vz]';Ҳ����û�м��ٶȲ���
    
    properties (SetAccess = protected)
        F;                  % ״̬ת�ƾ���
    end
    
    methods
        %% ���캯��
        % ���������
        % t�������������֮֡��ļ�࣬Ҳ����ʱ�䲽��
        % std�����������ı�׼��
        function td = CVmodel(t, std)
            td.T = t;
            td.Std = std;
            dimension = length(std);
            % ����״̬ת�ƾ���
            f = [1, t; 0, 1];
            o = zeros(2,2);
            td.F = [];
            for k = 1:dimension
                temp = [];
                for j = 1:dimension
                    if k == j
                        temp = [temp, f];
                    else
                        temp = [temp, o];
                    end
                end
                td.F = [td.F;temp];
            end
            % ���ɹ���������Э�������
            td.CovMatrix = [];
            q = [1/4 * t^4, 1/2 * t^3; 1/2 * t^3, t^2];
            for k = 1:dimension
                temp = [];
                for j = 1:dimension
                    if k == j
                        temp = [temp, q * std(k)^2];
                    else
                        temp = [temp, o];
                    end
                end
                td.CovMatrix = [td.CovMatrix;temp];
            end
        end
        %% һ�����ƺ���
        % ���һ��Ԥ��ֵ���ú���һ�����ڿ������˲���
        function x_new = estimateNewState(obj, x_old)
            x_new = obj.F * x_old;
        end
        % ��x������Ϸֲ��Ĺ���������һ�����������˲��е����ӳ�ʼ����
        function x_new = addDisturb(obj,x)
            b = [obj.T^2/2;obj.T];
            o = zeros(2,1);
            B = [];
            dimension = length(obj.Std);
            for k = 1:dimension
                temp = [];
                for j = 1:dimension
                    if k == j
                        temp = [temp, b];
                    else
                        temp = [temp, o];
                    end
                end
                B = [B;temp];
            end
            noise = B * diag(obj.Std) * randn(dimension,1);
            x_new = x + noise;
        end
        % ��Ԥ��ֵ�����������Ի�õ�ǰʱ��Ŀ��״̬������ֲ������������˲��С�
        % �ú���������ǰ���������ϳɵõ�
        function x_new = predictNewState(obj,x_old)
            b = [obj.T^2/2;obj.T];
            o = zeros(2,1);
            B = [];
            dimension = length(obj.Std);
            for k = 1:dimension
                temp = [];
                for j = 1:dimension
                    if k == j
                        temp = [temp, b];
                    else
                        temp = [temp, o];
                    end
                end
                B = [B;temp];
            end
            noise = B * diag(obj.Std) * randn(dimension,1);
            x_new = obj.F * x_old + noise;
        end
        %% ���״̬���̵����Խ���
        % ���������Ա任�����ֱ�����״̬ת�ƾ���
        function F_linear = getLinear(obj,~)
            F_linear = obj.F;
        end
    end
    
    methods (Static)
        % ��������λ�á��ٶȡ����ٶȵ���ʽ�ֱ����
        % �ֱ�Ҫ���Ƕ�ά����λ���
        % CVģ��û�м��ٶ���һ��������˸�ֵ��0
        function [pos,vel,acc] = transToCart(data, dimension)
            if dimension == 2
                pos = [data(1,:);data(3,:)];
                vel = [data(2,:);data(4,:)];
                acc = zeros(size(pos));
            elseif dimension == 3
                pos = [data(1,:);data(3,:);data(5,:)];
                vel = [data(2,:);data(4,:);data(6,:)];
                acc = zeros(size(pos));
            end
        end
    end
end

