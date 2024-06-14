classdef CVmodel < Model
    % 常速度模型
    % 状态量为：
    % state = [x, vx, y, vy, z, vz]';也就是没有加速度部分
    
    properties (SetAccess = protected)
        F;                  % 状态转移矩阵
    end
    
    methods
        %% 构造函数
        % 输入变量：
        % t：仿真过程中两帧之间的间距，也就是时间步长
        % std：过程噪声的标准差
        function td = CVmodel(t, std)
            td.T = t;
            td.Std = std;
            dimension = length(std);
            % 生成状态转移矩阵
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
            % 生成过程噪声的协方差矩阵
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
        %% 一步递推函数
        % 获得一步预测值，该函数一般用在卡尔曼滤波中
        function x_new = estimateNewState(obj, x_old)
            x_new = obj.F * x_old;
        end
        % 给x加入符合分布的过程噪声，一般用于粒子滤波中的粒子初始化中
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
        % 给预测值加噪声，可以获得当前时刻目标状态的先验分布，用于粒子滤波中。
        % 该函数可以由前两个函数合成得到
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
        %% 获得状态方程的线性近似
        % 由于是线性变换，因此直接输出状态转移矩阵
        function F_linear = getLinear(obj,~)
            F_linear = obj.F;
        end
    end
    
    methods (Static)
        % 将数据以位置、速度、加速度的形式分别输出
        % 分别要考虑二维和三位情况
        % CV模型没有加速度这一个量，因此该值赋0
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

