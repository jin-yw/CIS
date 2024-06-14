function [overlappingRate] = overlappingRateCalc(Range,Agent, TotalRangeX, AreaSize)
PointNum = 2000000;
PointCount = zeros(1,PointNum);
C = [];
for n = 1:size(TotalRangeX,1)
    C = [C; unifrnd(TotalRangeX(n,1), TotalRangeX(n,2),1,PointNum)];
end
AgentNum = length(Agent);
parfor idx = 1:PointNum
    flag = false;
    for i = 1:AgentNum
        pd = Agent{i}.getPd(C(:,idx));
        if pd > 0
            flag = true;
            break;
        end
    end
    if flag
        PointCount(idx) = 1;
    end
end
A_bar = AreaSize * sum(PointCount) / PointNum;
A = zeros(1,AgentNum);
parfor idx = 1:AgentNum
    tempR = Range{idx}(1,2);
    tempTheta = Range{idx}(2,2) - Range{idx}(2,1);
    A(idx) = pi * tempR^2 * tempTheta / (2*pi);
end
overlappingRate = (AgentNum * (sum(A) - A_bar) ) / ((AgentNum - 1) * sum(A));
end

