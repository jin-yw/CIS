classdef Cartesian < Model
    % 在直角坐标系下进行测量
    % 测量需要与状态相适应
    properties
        H;          % 测量矩阵
        StateType;  % 运动方程类型：CA模型、CT模型、CV模型
    end
    
    methods
        %% 构造函数
        function td = Cartesian(t,std,stateType)
            CVlike = {'CVmodel'};
            CAlike = {'CAmodel'};
            td.T = t;                   % 仿真周期：虽然测量方程中不会用到，但是基类中含有此变量，因此仍进行赋值
            td.Std = std;               % 测量噪声标准差，应该与测量的维数相适应
            td.StateType = stateType;
            dimension = length(std);
            % 计算得到测量矩阵
            switch stateType
                case CVlike
                    h = [1, 0];
                    o = zeros(1,2);
                case CAlike
                    h = [1, 0, 0];
                    o = zeros(1,3);
%                 case 'CTmodel'
%                     h = [1,0,0];
%                     o = zeros(1,3);
            end
            td.H = [];
            for k = 1:dimension
                temp = [];
                for j = 1:dimension
                    if k == j
                        temp = [temp, h];
                    else
                        temp = [temp, o];
                    end
                end
                td.H = [td.H;temp];
            end
            % 计算得到测量噪声的协方差矩阵
            td.CovMatrix = diag(std.*std);
        end
        %% 测量函数
        % 由于是在直角坐标系下测量，因此是一个线性的过程
        function  z_new = estimateNewState(obj, state)
            z_new = obj.H * state;
        end
        % 对于计算得到的测量值进行加噪
        function z_new = addDisturb(obj,z)
            dimension = length(obj.Std);
            z_new = z + chol(obj.CovMatrix) * randn(dimension, 1);
        end
        % 获得带噪声的测量值，由前两个函数合并实现
        function z_new = predictNewState(obj,state)
            temp = obj.estimateNewState(state);
            z_new = obj.addDisturb(temp);
        end
        %% 获得状态方程的线性近似
        % 由于是在直角坐标系下进行测量，且运动方程一般也是在直角坐标系中建立，不涉及直角坐标、球坐标之间相互转换的问题
        % 直接输出测量矩阵即可
        function H_linear = getLinear(obj, ~)
            H_linear = obj.H;
        end
        %% 获得概率
        % 该函数一般不会用在卡尔曼滤波中，在多目标滤波以及粒子滤波中会采用。
        % 用来获得p(y|x)
        function prob = getProb(obj, measure, state)
            z = obj.estimateNewState(state);
            e = z - measure;
%             if e' * e < obj.Std * obj.Std'
%                 e' * e
%             end
            prob = 1 / sqrt(det((2 * pi) * obj.CovMatrix)) * exp(-1/2 * e' / obj.CovMatrix * e);
        end
    end
    
end

