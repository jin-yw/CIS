classdef GeneratorUniform < handle
    % PHD�в��õ��Ӳ�ģ��
    
    properties
        Average;            % �Ӳ��ľ�ֵ��ʵ���ϣ������ʵ��ʱ����Ϊ0��ֻͨ�����÷ֲ���Χȷ��
        Range;              % �Ӳ��ֲ���Χ
    end
    
    methods
        %% ���캯��
        % ���������
        % average���Ӳ���ֵ����ʵ�ڱ��������ô��Ѿ���̫���ˣ�
        % range���Ӳ��ֲ���Χ
        function td = GeneratorUniform(average,range)
            td.Average = average;
            td.Range = range;
        end
        %% ����б���
        % �������б����չʾ��ĸ�������
        function returnList = giveList(obj)
            num = size(obj.Average,2);
            returnList = cell(num+1,1);
            for k = 1:num
                returnList{k} = ['̽���λ�ã�',mat2str(obj.Average(:,k))];
            end
            returnList{num+1} = ['̽�ⷶΧ��',mat2str(obj.Range)];
        end
        % ���ڷ�����ĸ������Ե��������ƣ����ڵ���ʽ�˵���
        function returnList = giveProperty(obj)
            returnList = {'̽���λ��','̽�ⷶΧ'};
        end
        % �����б���е�����λ�ã����ص���ʽ�˵��е���������λ��
        function returnIndex = givePropIndex(obj,index)
            num = size(obj.Average,2);
            if index <= num
                returnIndex = 1;
            else
                returnIndex = 2;
            end
        end
        % �����б���е�����λ�ã�������Ӧ���Ե�ֵ
        function returnValue = giveValue(obj,index)
            num = size(obj.Average,2);
            if index <= num
                returnValue = obj.Average(:,index);
            else
                returnValue = obj.Range;
            end
        end
        % ���뵯��ʽ�˵��е����ݣ������б������Ӧ���Ե�λ�ú�����ֵ������ֵ����ڲ���ֵ��Ӧ�Ŀɱ༭�ı�����
        function [returnValue,index] = giveValueStr(obj,str)
            switch str
                case '̽���λ��'
                    returnValue = obj.Average(:,1);
                    index = 1;
                case '̽��Э�������'
                    returnValue = obj.Range;
                    index = size(obj.Average,2)+1;
            end
        end
        %% Set����
        % ��������ֵ������indexʱ����λ���б���е�λ�ã�data�Ǹ��ĵ�����ֵ
        function setProperty(obj,index,data)
            switch index
                case 1
                    % ̽���λ��
                    obj.addAverage(data);
                case 2
                    % ̽��Э�������
                    obj.setMatrix(data);
            end
        end
        % ɾ������ֵ���ڱ�ģ���У�����ɾ��������ֵֻ�зֲ��ľ�ֵ
        % ��һ��˵��GeneratorUniform����Ըĳ�������GeneratorGaussian���Ӧ�ó���
        function delProperty(obj,index)
            num = size(obj.Average,2);
            if index <= num
                temp = zeros(size(obj.Average,1), num-1);
                for k = 1:num
                    if k < index
                        temp(:,k) = obj.Average(:,k);
                    elseif k > index
                        temp(:,k) = obj.Average(:,k-1);
                    end
                end
                obj.Average = temp;
            else
                return;
            end
        end
        % ���þ�ֵ��λ��
        function setAverage(obj,average)
%             if size(obj.Average,1) == size(average,1)
                obj.Average = average;
%             end
        end
        % ���Ӿ�ֵ��λ��
        function addAverage(obj,average)
            if size(obj.Average,1) == size(average,1)
                for k = 1:size(average,2)
                    flag = false;
                    for j = 1:size(obj.Average,2)
                        if isequal(obj.Average(:,j),average(:,k))
                            % ����ֻ���Ƿ�������ж�
                            flag = true;
                            break;
                        end
                    end
                    if ~flag
                        obj.Average = [obj.Average, average(:,j)];
                    end
                end
            end
        end
        % ���÷�Χ
        function setMatrix(obj,range)
            if size(range,1) == size(obj.Average, 1) && size(range,2) == size(obj.Range,2)
                obj.Range = range;
            end
        end
        %% ���ɵ�
        function point = generatePoint(obj)
            [~,num] = size(obj.Average);
            index = ceil(rand() * num);
            point = obj.Average(:,index) + unifrnd(obj.Range(:,1),obj.Range(:,2));
        end
        %% ��øõ�ĸ���
        function p = getProb(obj,~)
            range = obj.Range(:,2) - obj.Range(:,1);
            p = 1;
            for k = 1:length(range)
                if range(k) ~= 0
                    p = p * range(k);
                end
            end
            p = 1 / abs(p);
        end
    end
    
end

