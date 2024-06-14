classdef GeneratorUniform < handle
    % PHD中采用的杂波模型
    
    properties
        Average;            % 杂波的均值：实际上，在软件实现时设置为0，只通过设置分布范围确定
        Range;              % 杂波分布范围
    end
    
    methods
        %% 构造函数
        % 输入变量：
        % average：杂波均值（其实在本程序中用处已经不太大了）
        % range：杂波分布范围
        function td = GeneratorUniform(average,range)
            td.Average = average;
            td.Range = range;
        end
        %% 输出列表函数
        % 用于在列表框中展示类的各个参数
        function returnList = giveList(obj)
            num = size(obj.Average,2);
            returnList = cell(num+1,1);
            for k = 1:num
                returnList{k} = ['探测点位置：',mat2str(obj.Average(:,k))];
            end
            returnList{num+1} = ['探测范围：',mat2str(obj.Range)];
        end
        % 用于返回类的各个属性的中文名称，用于弹出式菜单中
        function returnList = giveProperty(obj)
            returnList = {'探测点位置','探测范围'};
        end
        % 输入列表框中的属性位置，返回弹出式菜单中的属性名称位置
        function returnIndex = givePropIndex(obj,index)
            num = size(obj.Average,2);
            if index <= num
                returnIndex = 1;
            else
                returnIndex = 2;
            end
        end
        % 输入列表框中的属性位置，返回相应属性的值
        function returnValue = giveValue(obj,index)
            num = size(obj.Average,2);
            if index <= num
                returnValue = obj.Average(:,index);
            else
                returnValue = obj.Range;
            end
        end
        % 输入弹出式菜单中的内容，返回列表框中相应属性的位置和属性值，属性值填充在参数值对应的可编辑文本框中
        function [returnValue,index] = giveValueStr(obj,str)
            switch str
                case '探测点位置'
                    returnValue = obj.Average(:,1);
                    index = 1;
                case '探测协方差矩阵'
                    returnValue = obj.Range;
                    index = size(obj.Average,2)+1;
            end
        end
        %% Set函数
        % 设置属性值，其中index时属性位于列表框中的位置，data是更改的属性值
        function setProperty(obj,index,data)
            switch index
                case 1
                    % 探测点位置
                    obj.addAverage(data);
                case 2
                    % 探测协方差矩阵
                    obj.setMatrix(data);
            end
        end
        % 删除属性值，在本模型中，可以删除的属性值只有分布的均值
        % 进一步说，GeneratorUniform类可以改成类似于GeneratorGaussian类的应用场景
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
        % 设置均值的位置
        function setAverage(obj,average)
%             if size(obj.Average,1) == size(average,1)
                obj.Average = average;
%             end
        end
        % 增加均值的位置
        function addAverage(obj,average)
            if size(obj.Average,1) == size(average,1)
                for k = 1:size(average,2)
                    flag = false;
                    for j = 1:size(obj.Average,2)
                        if isequal(obj.Average(:,j),average(:,k))
                            % 这里只用是否相等来判断
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
        % 设置范围
        function setMatrix(obj,range)
            if size(range,1) == size(obj.Average, 1) && size(range,2) == size(obj.Range,2)
                obj.Range = range;
            end
        end
        %% 生成点
        function point = generatePoint(obj)
            [~,num] = size(obj.Average);
            index = ceil(rand() * num);
            point = obj.Average(:,index) + unifrnd(obj.Range(:,1),obj.Range(:,2));
        end
        %% 获得该点的概率
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

