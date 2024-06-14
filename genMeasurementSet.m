function [measure] = genMeasurementSet(data, agent, Time, Lambda, range)
TimeLen = length(Time);
measure = cell(1, TimeLen);
trueState = cell(1, TimeLen);
dataLen = size(data,1);
for n = 1:dataLen
    tmpTraj = data{n, 2};
    tmpTrajLen = size(tmpTraj,1);
    for t = 1:tmpTrajLen
        currentTime = tmpTraj(t,1);
        index = find(Time == currentTime);
        x = [tmpTraj(t,2); 0; tmpTraj(t,3); 0];
        trueState{index} = [trueState{index}, x];
    end
end
%% Generate measurements
for k = 1:TimeLen
    TruthNum = size(trueState{k},2);
    if TruthNum > 0
        for n = 1:TruthNum
            x = trueState{k}(:,n);
            pd = agent.getPd(x);
            temp = rand();
            if temp <= pd
                measure{k} = [measure{k}, agent.MeasurementModel.predictNewState(x)];
            end
        end
        Nc = poissrnd(Lambda);
        C = [];
        for n = 1:size(range,1)
            C = [C; unifrnd(range(n,1), range(n,2), 1, Nc)];
        end
        measure{k} = [measure{k}, C];
    else
        Nc = poissrnd(Lambda);
        C = [];
        for n = 1:size(range,1)
            C = [C; unifrnd(range(n,1), range(n,2), 1, Nc)];
        end
        measure{k} = [measure{k}, C];
    end
end
end

