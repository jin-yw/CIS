classdef Cartesian < Model
    % ��ֱ������ϵ�½��в���
    % ������Ҫ��״̬����Ӧ
    properties
        H;          % ��������
        StateType;  % �˶��������ͣ�CAģ�͡�CTģ�͡�CVģ��
    end
    
    methods
        %% ���캯��
        function td = Cartesian(t,std,stateType)
            CVlike = {'CVmodel'};
            CAlike = {'CAmodel'};
            td.T = t;                   % �������ڣ���Ȼ���������в����õ������ǻ����к��д˱���������Խ��и�ֵ
            td.Std = std;               % ����������׼�Ӧ���������ά������Ӧ
            td.StateType = stateType;
            dimension = length(std);
            % ����õ���������
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
            % ����õ�����������Э�������
            td.CovMatrix = diag(std.*std);
        end
        %% ��������
        % ��������ֱ������ϵ�²����������һ�����ԵĹ���
        function  z_new = estimateNewState(obj, state)
            z_new = obj.H * state;
        end
        % ���ڼ���õ��Ĳ���ֵ���м���
        function z_new = addDisturb(obj,z)
            dimension = length(obj.Std);
            z_new = z + chol(obj.CovMatrix) * randn(dimension, 1);
        end
        % ��ô������Ĳ���ֵ����ǰ���������ϲ�ʵ��
        function z_new = predictNewState(obj,state)
            temp = obj.estimateNewState(state);
            z_new = obj.addDisturb(temp);
        end
        %% ���״̬���̵����Խ���
        % ��������ֱ������ϵ�½��в��������˶�����һ��Ҳ����ֱ������ϵ�н��������漰ֱ�����ꡢ������֮���໥ת��������
        % ֱ������������󼴿�
        function H_linear = getLinear(obj, ~)
            H_linear = obj.H;
        end
        %% ��ø���
        % �ú���һ�㲻�����ڿ������˲��У��ڶ�Ŀ���˲��Լ������˲��л���á�
        % �������p(y|x)
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

