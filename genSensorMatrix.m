function [MeasurementModel, SensorPos, SensorNum, startTime, endTime, adjacentMatrix, Range] = genSensorMatrix(data, radius, std, Step)
%% Find the minimum and maximum of the position
dataLen = size(data,1);
minX = min(data{1,2}(:,2));
minX = minX(1);
maxX = max(data{1,2}(:,2));
maxX = maxX(1);
minY = min(data{1,2}(:,3));
minY = minY(1);
maxY = max(data{1,2}(:,3));
maxY = maxY(1);
minT = min(data{1,2}(:,1));
minT = minT(1);
maxT = max(data{1,2}(:,1));
maxT = maxT(1);
for idx = 2:dataLen
    tempMinX = min(data{idx,2}(:,2));
    tempMaxX = max(data{idx,2}(:,2));
    tempMinY = min(data{idx,2}(:,3));
    tempMaxY = max(data{idx,2}(:,3));
    tempMinT = min(data{idx,2}(:,1));
    tempMaxT = max(data{idx,2}(:,1));
    tempMinX = tempMinX(1);
    tempMaxX = tempMaxX(1);
    tempMinY = tempMinY(1);
    tempMaxY = tempMaxY(1);
    tempMinT = tempMinT(1);
    tempMaxT = tempMaxT(1);
    if tempMinX < minX
        minX = tempMinX;
    end
    if tempMaxX > maxX
        maxX = tempMaxX;
    end
    if tempMinY < minY
        minY = tempMinY;
    end
    if tempMaxY > maxY
        maxY = tempMaxY;
    end
    if tempMinT < minT
        minT = tempMinT;
    end
    if tempMaxT > maxT
        maxT = tempMaxT;
    end
end
startTime = minT;
endTime = maxT;
dX = maxX - minX;
dY = maxY - minY;
Xnum = ceil(round(dX) / radius) + 1;
Ynum = ceil(round(dY) / radius) + 1;
SensorNum = Xnum * Ynum;
SensorPos = cell(1, SensorNum);
adjacentMatrix = zeros(SensorNum, SensorNum);
for idx = 1:SensorNum
    row = ceil(idx / Xnum);
    col = idx - (row-1)*Xnum;
    posX = minX + (col-1) * radius;
    posY = maxY - (row-1) * radius;
    SensorPos{idx} = [posX; posY];
    adjacentMatrix(idx, idx) = 1;
    if col == 1
        if row == 1
            adjacentMatrix(idx, idx + 1) = 1;
            adjacentMatrix(idx, idx + Xnum) = 1;
        elseif row == Ynum
            adjacentMatrix(idx, idx + 1) = 1;
            adjacentMatrix(idx, idx - Xnum) = 1;
        else
            adjacentMatrix(idx, idx + 1) = 1;
            adjacentMatrix(idx, idx + Xnum) = 1;
            adjacentMatrix(idx, idx - Xnum) = 1;
        end
    elseif col == Xnum
        if row == 1
            adjacentMatrix(idx, idx - 1) = 1;
            adjacentMatrix(idx, idx + Xnum) = 1;
        elseif row == Ynum
            adjacentMatrix(idx, idx - 1) = 1;
            adjacentMatrix(idx, idx - Xnum) = 1;
        else
            adjacentMatrix(idx, idx - 1) = 1;
            adjacentMatrix(idx, idx - Xnum) = 1;
            adjacentMatrix(idx, idx + Xnum) = 1;
        end
    else
        if row == 1
            adjacentMatrix(idx, idx - 1) = 1;
            adjacentMatrix(idx, idx + 1) = 1;
            adjacentMatrix(idx, idx + Xnum) = 1;
        elseif row == Ynum
            adjacentMatrix(idx, idx - 1) = 1;
            adjacentMatrix(idx, idx + 1) = 1;
            adjacentMatrix(idx, idx - Xnum) = 1;
        else
            adjacentMatrix(idx, idx - 1) = 1;
            adjacentMatrix(idx, idx + 1) = 1;
            adjacentMatrix(idx, idx - Xnum) = 1;
            adjacentMatrix(idx, idx + Xnum) = 1;
        end
    end
end
MeasurementModel = cell(1, SensorNum);
Range = cell(1, SensorNum);
for idx = 1:SensorNum
    MeasurementModel{idx} = SphereSensor(Step, std, SensorPos{idx}, 'CVmodel');
    Range{idx} = [0, radius; -pi, pi];
end
end

