classdef GeneratorGaussian < handle
    % PHD�в��õ�Ŀ������
    
    properties
        Average;            % ���˹�ֲ���ÿ����˹�ֲ���λ��
        Covarience;         % Э�������
        Prob;               % The intensity of the distribution.
    end
    
    methods
        %% ���캯��
        % ���������
        % average�����˹�ֲ���ÿ����˹�ֲ���λ��
        % range��Э�������
        function td = GeneratorGaussian(average, range, prob)
            if nargin == 2
                td.Average = average;
                td.Covarience = range;
                td.Prob = 0.2;
            elseif nargin == 3
                td.Average = average;
                td.Covarience = range;
                td.Prob = prob;
            end
        end
        %% Generate a point
        function point = generatePoint(obj)
            [dim,num] = size(obj.Average);
            index = ceil(rand() * num);
            matrix = chol(obj.Covarience,'lower');
            point = obj.Average(:,index) + matrix * randn(dim,1);
        end
        %% Get the probability of the point
        function p = getProb(obj,point)
            [dim,num] = size(obj.Average);
            p = 0;
            for k = 1:num
                tempP = 1 / (sqrt(2 * pi)^(dim) * det(obj.Covarience)^(1/2) ) * exp(-1/2 * (obj.Average(:,k) - point)' / obj.Covarience * (obj.Average(:,k) - point));
                p = p + tempP / num;
            end
        end
        %% Output the 1st and 2nd moment of Gaussian distribution
        function model = generateModel(obj)
            model.x = obj.Average;
            model.P = obj.Covarience;
            model.w = obj.Prob;
        end
    end
    
end

