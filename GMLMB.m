classdef GMLMB < FilterAlgorithm
    % The LMB filter with adaptive birth distribution.
    
    properties
        BernoulliSet;
        GenerateSet;
        ClutterModel;
        Pg;
        Ps;
        Time;
        Gamma;
        Lambda;
        LambdaB;
        rBMax;
        THETA;
    end
    
    methods
        % Constructor
        function td = GMLMB(kinematicmodel, agent, cluttermodel, ps, pg, gamma, lambda, lambdaB, rbmax, theta)
            td.BernoulliSet = cell(0,0);
            td.GenerateSet = cell(0,0);
            td.KinematicModel = kinematicmodel;
            td.Agent = agent;
            td.ClutterModel = cluttermodel;
            td.Pg = pg;
            td.Ps = ps;
            td.Time = 1;
            td.Gamma = gamma;
            td.Lambda = lambda;
            td.LambdaB = lambdaB;
            td.rBMax = rbmax;
            td.THETA = theta;
        end
        % The batch filter
        function [Track, Distribution] = filter(this, measure)
            timeLength = length(measure);
            Track = cell(0,0);
            this.Time = 1;
            this.BernoulliSet = cell(0,0);
            this.GenerateSet = cell(0,0);
            Distribution = cell(1,timeLength);
            for k = 1:timeLength
                time = k;
%                 time
                [tempBernoulliSet, Distribution{k}] = this.filter_update(measure{k});
                Track = this.updateTrack(tempBernoulliSet, Track);
                this.Time = this.Time + 1;
            end
        end
        function init(this)
            this.Time = 1;
            this.BernoulliSet = cell(0,0);
            this.GenerateSet = cell(0,0);
        end
        function TimeUpdate(this)
            this.Time = this.Time + 1;
        end
        function Track = genTrack(this, LMB)
            timeLength = length(LMB);
            Track = cell(0,0);
            this.Time = 1;
            for k = 1:timeLength
                Track = this.updateTrack(LMB{k}, Track);
                this.Time = this.Time + 1;
            end
        end
        % The filter for single step
        function [outputBernoulliSet, BernoulliDistribution] = filter_update(this, measure)
            %% Prediction
            oldBernoulliLen = length(this.BernoulliSet);
            newbornLen = length(this.GenerateSet);
            preBernoulliSetOld = cell(1,oldBernoulliLen);
            preBernoulliSetNewborn = cell(1,newbornLen);
            % prediction for existed targets
            parfor idx = 1:oldBernoulliLen
                preBernoulliSetOld{idx} = this.updateBernoulliComponent(this.BernoulliSet{idx});
            end
            % prediction for newborn targets
            parfor idx = 1:newbornLen
                preBernoulliSetNewborn{idx} = this.genNewbornTarget(idx);
            end
            preBernoulliSet = this.combineBernoulliSet(preBernoulliSetOld, preBernoulliSetNewborn);
            %% Grouping and Gating
            LMBgroup = this.genGroup(preBernoulliSet, measure);
            %% Parallel Group Update
            % The predicted LMBs are represented as the form of \delta-GLMB
            LMBgroupLen = length(LMBgroup);
            preGLMB = cell(1,LMBgroupLen);
            parfor idx = 1:LMBgroupLen
                preGLMB{idx} = this.LMB2GLMB(LMBgroup{idx}, 30);
            end
            % Parallel \delta-GLMB Update
            estGLMB = cell(1,LMBgroupLen);
            estLMBset = cell(1,LMBgroupLen);
            parfor idx = 1:LMBgroupLen
                estGLMB{idx} = this.getUpdateGLMB(preGLMB{idx}, measure(:, LMBgroup{idx}.Z));
                estLMBset{idx} = this.GLMB2LMB(estGLMB{idx});
            end
            preGenerateSet = this.genBirthDistribution(estGLMB, LMBgroup, measure);
            % Assemble LMB
            LMB = this.assembleLMB(estLMBset);
            %% Merging
            LMB = this.trackMerging(LMB);
            %% Track pruning
            LMB = this.trackPruning(LMB);
            this.GenerateSet = this.trackPruning(preGenerateSet);
            %% Target estimation
            N = this.getN(LMB);
%             N
            outputBernoulliSet = cell(1,N);
            parfor idx = 1:N
                outputBernoulliSet{idx} = LMB{idx};
            end
            this.BernoulliSet = LMB;
            BernoulliDistribution = LMB;
        end
    end
    
    methods
        function newTrack = updateTrack(this, BernoulliSet, Track)
            newTrack = Track;
            TrackLen = length(Track);
            BernoulliSetLen = length(BernoulliSet);
            count = TrackLen + 1;
            for idx1 = 1:BernoulliSetLen
                flag = true;
                for idx2 = 1:TrackLen
                    if isequal(newTrack{idx2}.l, BernoulliSet{idx1}.l)
                        tempX = BernoulliSet{idx1}.x;
                        newTrack{idx2}.X = [newTrack{idx2}.X, tempX];
                        newTrack{idx2}.t = [newTrack{idx2}.t, this.Time];
                        flag = false;
                        break;
                    end
                end
                if flag
                    tempX = BernoulliSet{idx1}.x;
                    newTrack{count}.X = tempX;
                    newTrack{count}.t = this.Time;
                    newTrack{count}.l = BernoulliSet{idx1}.l;
                    count = count + 1;
                end
            end
        end
        % prediction for existed targets:
        function bernoulliComponent = updateBernoulliComponent(this, oldBernoulliComponent)
            bernoulliComponent.r = this.Ps * oldBernoulliComponent.r;
            bernoulliComponent.l = oldBernoulliComponent.l;
            bernoulliComponent.x = this.KinematicModel.estimateNewState(oldBernoulliComponent.x);
            F = this.KinematicModel.getLinear(oldBernoulliComponent.x);
            bernoulliComponent.P = F * oldBernoulliComponent.P * F' + this.KinematicModel.CovMatrix;
        end
        % prediction for newborn targets
        function bernoulliComponent = genNewbornTarget(this, index)
            bernoulliComponent.x = this.GenerateSet{index}.x;
            bernoulliComponent.P = this.GenerateSet{index}.P;
            bernoulliComponent.r = this.GenerateSet{index}.r;
            bernoulliComponent.l = [this.Time, index];
        end
        % combine the Bernoulli sets
        function BernoulliSet = combineBernoulliSet(this, BernoulliSet1, BernoulliSet2)
            BernoulliSet1Len = length(BernoulliSet1);
            BernoulliSet2Len = length(BernoulliSet2);
            BernoulliSet = cell(1,BernoulliSet1Len + BernoulliSet2Len);
            parfor idx = 1:BernoulliSet1Len
                BernoulliSet{idx} = BernoulliSet1{idx};
            end
            parfor idx = 1:BernoulliSet2Len
                BernoulliSet{idx + BernoulliSet1Len} = BernoulliSet2{idx};
            end
        end
        % grouping and gating
        function LMBset = genGroup(this, BernoulliSet, measure)
            % Initialization
            BernoulliSetLen = length(BernoulliSet);
            LMBset = cell(1,BernoulliSetLen);
            measureLen = size(measure,2);
            for i = 1:BernoulliSetLen
                tempLMB.BernoulliComponent = cell(0,0);
                tempLMB.BernoulliComponent{1} = BernoulliSet{i};
                tempLMB.Z = [];
                for j = 1:measureLen
                    flag = this.ifInRange(BernoulliSet{i}.x, BernoulliSet{i}.P, measure(:,j));
                    if flag
                        tempLMB.Z = [tempLMB.Z, j];
                    end
                end
                LMBset{i} = tempLMB;
            end
            % Iteration
            LMBset = this.genGroupIteration(LMBset);
        end
        function flag = ifInRange(this, x, P, z)
            H = this.Agent.MeasurementModel.getLinear(x);
            S = H * P * H' + this.Agent.MeasurementModel.CovMatrix;
            zz = this.Agent.MeasurementModel.estimateNewState(x);
            temp = (z - zz)' / S * (z - zz);
%             temp = sqrt(temp);
            if temp < this.Gamma
                flag = true;
            else
                flag = false;
            end
        end
        function [newLMBset] = genGroupIteration(this, oldLMBset)
            oldLMBsetLen = length(oldLMBset);
            newLMBset = cell(0,0);
            for n = 1:oldLMBsetLen
                newLMBset = this.appendLMBset(newLMBset, oldLMBset{n});
            end
            if length(newLMBset) < length(oldLMBset)
                newLMBset = this.genGroupIteration(newLMBset);
            end
        end
        function [newLMBset] = appendLMBset(this, oldLMBset, LMBset)
            oldLMBsetLen = length(oldLMBset);
            newLMBset = oldLMBset;
            flag = false;
            for idx = 1:oldLMBsetLen
                tempClass = intersect(oldLMBset{idx}.Z, LMBset.Z);
                if ~isempty(tempClass)
                    newLMBset{idx}.BernoulliComponent = [newLMBset{idx}.BernoulliComponent, LMBset.BernoulliComponent];
                    newLMBset{idx}.Z = union(newLMBset{idx}.Z, LMBset.Z);
                    flag = true;
                    break;
                end
            end
            if ~flag
                newLMBset{oldLMBsetLen + 1} = LMBset;
            end
        end
        % LMB RFS to GLMB RFS
        function [GLMB] = LMB2GLMB(this, LMBset, Num)
            bernoulliSet = LMBset.BernoulliComponent;
            bernoulliSetLen = length(bernoulliSet);
            r = zeros(bernoulliSetLen,1);
            Distribution = cell(1,bernoulliSetLen);
            Label = zeros(bernoulliSetLen,2);
            parfor idx = 1:bernoulliSetLen
                Distribution{idx}.x = bernoulliSet{idx}.x;
                Distribution{idx}.P = bernoulliSet{idx}.P;
                Label(idx,:) = bernoulliSet{idx}.l;
                r(idx) = bernoulliSet{idx}.r;
            end
            GLMB.Distribution = Distribution;
            GLMB.l = Label;
            cost = r./(1-r);
            neglogcost = -log(cost);
            [paths,nlcost] = this.ksp(neglogcost, Num);
            nlcostLen = length(nlcost);
            GLMB.w = zeros(1,nlcostLen);
            GLMB.I = cell(1,nlcostLen);
            for idx = 1:nlcostLen
                GLMB.w(idx) = sum(log(1-r)) - nlcost(idx);
                GLMB.I{idx} = paths{idx};
            end
            GLMB.w = exp(GLMB.w);
            GLMB.w = GLMB.w / sum(GLMB.w);
        end
        function [paths, nlcost] = ksp(this, neglogcost, Num)
            if Num == 0
                paths = [];
                nlcost = [];
            else
                % P = empty
                % P: set of shortest paths from s to u
                P = cell(0,0);
                % count_u = 0 for all u in V
                % count_u: number of shortest paths found to node u
                % V: set of vertices
                [graph, minValue] = this.genGraph1(neglogcost);
                graphLen = length(graph);
                count_u = cell(1,graphLen);
                for idx = 1:graphLen
                    tempLen = length(graph{idx});
                    count_u{idx} = zeros(tempLen,1);
                end
                % insert path Ps={s} into B with cost 0
                P_s = [1;1];
                B{1} = P_s;
                cost = graph{1}(1);
                nlcost = [];
                while ~isempty(B) && length(P) < Num
                    % let Pu be the shortest cost path in B with cost C
                    minCost = min(cost);
                    index = find(cost == minCost);
                    index = index(1);
                    Pu = B{index};
                    C = minCost;
                    u = Pu(:,end);
                    % B = B-{Pu}, count_u = count_u + 1;
                    newB = cell(1,length(cost)-1);
                    newCost = zeros(1,length(cost)-1);
                    for idx = 1:length(cost)
                        if idx < index
                            newB{idx} = B{idx};
                            newCost(idx) = cost(idx);
                        elseif idx > index
                            newB{idx-1} = B{idx};
                            newCost(idx-1) = cost(idx);
                        end
                    end
                    B = newB;
                    cost = newCost;
                    count_u{u(1)}(u(2)) = count_u{u(1)}(u(2)) + 1;
                    % if u = t then P = P \cup {Pu}
                    if u(1) == graphLen
                        PLen = length(P);
                        P{PLen + 1} = Pu;
                        nlcost = [nlcost,C];
                        continue;
                    end
                    % if count_u <= K then
                    if count_u{u(1)}(u(2)) <= Num
                        % for each vertex v adjacent u:
                        % let Pv be a new path with cost C+w(u,v) formed by
                        % concatenating edge (u,v) to path Pu
                        nextLayer = u(1) + 1;
                        tempLen = length(graph{nextLayer});
                        Pv = cell(1,tempLen);
                        costv = zeros(1,tempLen);
                        for idx = 1:tempLen
                            Pv{idx} = [Pu, [nextLayer; idx]];
                            costv(idx) = C + graph{nextLayer}(idx);
                        end
                        % insert Pv into B
                        BLen = length(cost);
                        for idx = 1:tempLen
                            B{BLen + idx} = Pv{idx};
                            cost(BLen + idx) = costv(idx);
                        end
                    end
                end
                PLen = length(P);
                paths = cell(1,PLen);
                for idx = 1:PLen
                    tempPath = P{idx};
                    for i = 1:graphLen
                        if i == 1 || i == graphLen
                            continue;
                        else
                            tempIndex = tempPath(2,i);
                            if tempIndex == 1
                                paths{idx} = [paths{idx}, i-1];
                            end
                        end
                    end
                end
                nlcost = nlcost + minValue * graphLen;
            end
        end
        function [graph, minValue] = genGraph1(this, value)
            valueLen = length(value);
            Len = valueLen + 2;
            graph = cell(1,Len);
            minValue = 0;
            for idx = 1:Len
                if idx == 1
                    graph{idx} = 0;
                elseif idx == Len
                    graph{idx} = 0;
                else
                    graph{idx} = [value(idx-1); 0];
                    if minValue > value(idx-1)
                        minValue = value(idx-1);
                    end
                end
            end
            parfor idx1 = 1:Len
                graph{idx1} = graph{idx1} - minValue;
            end
        end
        % Get the updated GLMB distribution
        function [GLMB] = getUpdateGLMB(this, preGLMB, measure)
            measureLen = size(measure,2);
            % update each Benroulli Component
            BernoulliLen = length(preGLMB.Distribution);
            % Correlated with no measurements
            [Distribution, Label, ETA] = this.updateBernoulli2(preGLMB, measure, 0, BernoulliLen);
            Z = zeros(1,BernoulliLen);
            % Correlated with idx_th measurement
            for idx = 1:measureLen
                [tempDistribution, tempLabel, tempETA] = this.updateBernoulli2(preGLMB, measure, idx, BernoulliLen);
                Distribution = [Distribution, tempDistribution];
                Label = [Label; tempLabel];
                ETA = [ETA; tempETA];
                Z = [Z; ones(1,BernoulliLen) * idx];
            end
            Z = Z';
            % Update the GLMB weight
            I_len = length(preGLMB.I);
            posteriorGLMB = cell(1,I_len);
            parfor idx = 1:I_len
                posteriorGLMB{idx} = this.getPosteriorGLMB(preGLMB.l(preGLMB.I{idx},:), preGLMB.w(idx), Distribution, Label, ETA, Z, 20);
            end
            GLMB = this.assembleGLMB(posteriorGLMB, Distribution, Label, Z);
            GLMB.w = GLMB.w / sum(GLMB.w);
        end
        function [Distribution, Label, ETA] = updateBernoulli2(this, GLMB, measure, index, BernoulliLen)
            Distribution = cell(1,BernoulliLen);
            Label = zeros(BernoulliLen,2);
            ETA = zeros(1, BernoulliLen);
            parfor idx2 = 1:BernoulliLen
                [Distribution{idx2}, eta] = this.updateBernoulli1(GLMB.Distribution{idx2}, measure, index);
                Label(idx2,:) = GLMB.l(idx2,:);
                ETA(idx2) = eta;
            end
        end
        function [newDistribution, eta] = updateBernoulli1(this, oldDistribution, measure, measureIdx)
            if measureIdx == 0
                eta = 1 - this.Agent.getPd(oldDistribution.x) * this.Pg;
                newDistribution = oldDistribution;
            else
                P_pre = oldDistribution.P;
                x_pre = oldDistribution.x;
                H = this.Agent.MeasurementModel.getLinear(x_pre);
                R = this.Agent.MeasurementModel.CovMatrix;
                S = H*P_pre*H' + R;
                K = P_pre * H' / S;
                z = this.Agent.MeasurementModel.estimateNewState(x_pre);
                dz = measure(:,measureIdx) - z;
                x_est = x_pre + K * dz;
                P_est = P_pre - K * H * P_pre;
                newDistribution.x = x_est;
                newDistribution.P = P_est;
                temp = 1/sqrt(det(2*pi*S)) * exp(-1/2 * dz' / S * dz);
                eta = this.Agent.getPd(oldDistribution.x) * this.Pg * temp / (this.ClutterModel.getProb(z) * this.Lambda);
            end
        end
        function newCost = correctCost(this, path, cost, minValue)
            costLen = length(cost);
            newCost = zeros(1, costLen);
            parfor j = 1:costLen
                tempPath = path{j};
                tempPathLen = size(tempPath,2);
                value = 0;
                for k = 1:tempPathLen
                    value = value + minValue;
                end
                newCost(j) = cost(j) + value;
            end
        end
        function [GLMB] = getPosteriorGLMB(this, LabelSet, w, Distribution, Label, ETA, Z, Num)
            % LabelSet: The label set for a specific GLMB component
            GLMB.Distribution = Distribution;
            GLMB.Label = Label;
            LabelSetLen = size(LabelSet, 1);
            if LabelSetLen == 0
                % No target assumption
                GLMB.w = w;
                GLMB.I = {[]};
            else
                indexSet = zeros(1,LabelSetLen);
                parfor idx1 = 1:LabelSetLen
                    for idx2 = 1:size(Label,1)
                        if isequal(LabelSet(idx1,:), Label(idx2,:))
                            indexSet(idx1) = idx2;
                            break;
                        end
                    end
                end
                tempETA = ETA(:,indexSet);
                C1 = tempETA(2:end,:);
                C1 = C1';
                [measureNum, targetNum] = size(ETA);
                measureNum = measureNum - 1;
                C1 = -log(C1);
                C2 = diag(tempETA(1,:));
                C2 = -log(C2);
                C = [C1, C2];
                minValue = min(min(C));
                if minValue < 0
                    C = C - minValue;
                end
                [cost, path] = MurtyAlgorithm(C, Num);
                costLen = length(cost);
                newPath = cell(1,costLen);
                parfor j = 1:costLen
                    [idx1, idx2] = find(path{j} == 1);
                    tempLen = length(idx1);
                    P = [];
                    for k = 1:tempLen
                        if idx2(k) > measureNum
                            node = [idx1(k);0];
                        else
                            node = [idx1(k); idx2(k)];
                        end
                        P = [P, node];
                    end
                    newPath{j} = P;
                end
                path = newPath;
                if minValue < 0
                    cost = this.correctCost(path, cost, minValue);
                end
                pathLen = length(path);
                w_est = zeros(pathLen,1);
                I = cell(1,pathLen);
                parfor i = 1:pathLen
                    tempPath = path{i};
                    tempPathLen = size(tempPath,2);
                    tempI = [];
                    w_est(i) = w * exp(-cost(i));
                    for j = 1:tempPathLen
                        if tempPath(2,j) ~= 0
                            I_index = targetNum * tempPath(2,j) + indexSet(tempPath(1,j));
                        else
                            I_index = indexSet(tempPath(1,j));
                        end
                        tempI = [tempI, I_index];
                    end
                    I{i} = tempI;
                end
                GLMB.I = I;
                GLMB.w = w_est;
            end
        end
        % Get the index of the Bernoulli component given label
        function index = getBerIndex1(this, label, BernoulliSet)
            BerLen = length(BernoulliSet);
            indexIndicator = zeros(1,BerLen);
            parfor idx = 1:BerLen
                if isequal(label, BernoulliSet{idx}.l)
                    indexIndicator(idx) = 1;
                end
            end
            index = find(indexIndicator == 1);
        end
        function index = getBerIndex2(this, label, LabelSet)
            Len = size(LabelSet, 1);
            indexIndicator = zeros(1,Len);
            parfor idx = 1:Len
                if isequal(label, LabelSet(idx,:))
                    indexIndicator(idx) = 1;
                end
            end
            index = find(indexIndicator == 1);
        end
        function indexSet = getAllCombination(this, label2dist, Z)
            label2dist_len = length(label2dist);
            tempSet1 = cell(0,0);
            for idx = 1:label2dist_len
                if idx == 1
                    Len1 = length(label2dist{idx});
                    for idx2 = 1:Len1
                        tempSet1{idx2} = label2dist{idx}(idx2);
                    end
                else
                    newTempSet = cell(0,0);
                    count = 1;
                    tempSet1Len = length(tempSet1);
                    for idx2 = 1:tempSet1Len
                        z1 = Z(tempSet1{idx2});
                        Z1 = Z(label2dist{idx});
                        for idx3 = 1:length(Z1)
                            if ~isempty(find(z1 == Z1(idx3))) && Z1(idx3) ~= 0
                                continue;
                            else
                                newTempSet{count} = [tempSet1{idx2}, label2dist{idx}(idx3)];
                                count = count + 1;
                            end
                        end
                    end
                    tempSet1 = newTempSet;
                end
            end
            indexSet = tempSet1;
        end
        function [Distribution, Label, ETA, Z] = pruneByETA(this, dist, l, Eta, BernoulliLen, measureLen, z)
            MaxNum = 40;
            minPercentage = 0.00001;
            if BernoulliLen * (measureLen + 1) <= MaxNum
                Distribution = dist;
                Label = l;
                ETA = Eta;
                Z = z;
            else
                MaxEta = max(Eta);
                minValue = MaxEta * minPercentage;
                tempIndex = find(Eta >= minValue);
                MaxNum = length(tempIndex);
                [~, index] = sort(Eta, 'descend');
                index = index(1:MaxNum);
                ETA = Eta(index);
                Distribution = cell(1,length(index));
                parfor idx = 1:length(index)
                    Distribution{idx} = dist{index(idx)};
                end
                Label = l(index,:);
                Z = z(index);
            end
        end
        % Re-assemble the posterior GLMB distributions
        function GLMB = assembleGLMB(this, oldGLMB, Distribution, Label, Z)
            GLMB.Distribution = Distribution;
            GLMB.Label = Label;
            GLMB.Z = Z;
            GLMB.I = {};
            GLMB.w = [];
            Len = length(oldGLMB);
            for idx = 1:Len
                if ~isempty(oldGLMB{idx})
                    GLMB.I = [GLMB.I, oldGLMB{idx}.I];
                    GLMB.w = [GLMB.w; oldGLMB{idx}.w];
                end
            end
        end
        % the \delta-GLMB form for each group is collapsed back to LMB
        % form.
        function LMB = GLMB2LMB(this, GLMB)
            Len = length(GLMB.Distribution);
            LMB = cell(1,Len);
            parfor idx = 1:Len
                LMB{idx}.Distribution = GLMB.Distribution{idx};
                LMB{idx}.r = 0;
                LMB{idx}.l = GLMB.Label(idx,:);
                LMB{idx}.z = GLMB.Z(idx);
            end
            wLen = length(GLMB.w);
            for idx = 1:wLen
                I = GLMB.I{idx};
                for idx1 = 1:length(I)
                    LMB{I(idx1)}.r = LMB{I(idx1)}.r + GLMB.w(idx);
                end
            end
        end
        % Obtain the multi-Benroulli birth distribution at the next time
        % step k+1 depends on the set of measurements Zk of the current
        % time step.
        function [preGenerateSet] = genBirthDistribution(this, GLMBset, LMBgroup, measure)
            measureLen = size(measure,2);
            LMBgroupLen = length(LMBgroup);
            rU = zeros(1,measureLen);
            for idx = 1:LMBgroupLen
                Z = LMBgroup{idx}.Z;
                GLMB = GLMBset{idx};
                ZLen = length(Z);
                for z = 1:ZLen
                    temp = 0;
                    I = GLMB.I;
                    ILen = length(I);
                    for i = 1:ILen
                        tempZ = GLMB.Z(I{i});
                        index = find(tempZ == z);
                        if ~isempty(index)
                            temp = temp + GLMB.w(i);
                        end
                    end
                    rU(Z(z)) = temp;
                end
            end
            rB = 1-rU;
            rB = rB / sum(rB);
            preGenerateSet = cell(1, measureLen);
            parfor idx = 1:measureLen
                preGenerateSet{idx}.x = this.genSingleTargetNewbornState(measure(:,idx));
                preGenerateSet{idx}.P = diag([50; 50; 50; 50].^2);
                preGenerateSet{idx}.r = min(this.rBMax, rB(idx) * this.LambdaB);
            end
        end
        function x = genSingleTargetNewbornState(this, z)
            pos = this.Agent.MeasurementModel.getCartesianPos(z);
            dimension = length(pos);
            CVlike = {'CVmodel'};
            CAlike = {'CAmodel'};
            switch class(this.KinematicModel)
                case CVlike
                    x = [];
                    for idx = 1:dimension
                        tempX = [pos(idx);0];
                        x = [x; tempX];
                    end
                case CAlike
                    x = [];
                    for idx = 1:dimension
                        tempX = [pos(idx);zeros(2,1)];
                        x = [x; tempX];
                    end
            end
        end
        % Re-assemble the posterior LMB distributions
        function LMB = assembleLMB(this, oldLMB)
            Len = length(oldLMB);
            LMB = {};
            for idx = 1:Len
                LMB = [LMB, oldLMB{idx}];
            end
        end
        % Merging
        function [newLMB] = trackMerging(this, LMB)
            LMBLen = length(LMB);
            newLMB = cell(0,0);
            count = 1;
            for idx = 1:LMBLen
                index = this.getBerIndex1(LMB{idx}.l, newLMB);
                if isempty(index)
                    if LMB{idx}.r > 0
                        newLMB{count}.l = LMB{idx}.l;
                        newLMB{count}.x = LMB{idx}.Distribution.x;
                        newLMB{count}.P = LMB{idx}.Distribution.P;
                        newLMB{count}.r = LMB{idx}.r;
                        count = count + 1;
                    end
                else
                    index = index(1);
                    x1 = newLMB{index}.x;
                    x2 = LMB{idx}.Distribution.x;
                    P1 = newLMB{index}.P;
                    P2 = LMB{idx}.Distribution.P;
                    r1 = newLMB{index}.r;
                    r2 = LMB{idx}.r;
                    P = inv((r1*inv(P1) + r2*inv(P2))/(r1+r2));
                    x = P * (r1*inv(P1)*x1 + r2*inv(P2)*x2)/(r1+r2);
                    r = r1 + r2;
                    newLMB{index}.x = x;
                    newLMB{index}.P = P;
                    newLMB{index}.r = r;
                end
            end
            newLMBLen = length(newLMB);
            oldLMB = newLMB;
            newLMB = cell(1,newLMBLen);
            parfor idx = 1:newLMBLen
                if oldLMB{idx}.r > 1
                    newLMB{idx} = oldLMB{idx};
                    newLMB{idx}.r = 1;
                else
                    newLMB{idx} = oldLMB{idx};
                end
            end
        end
        % Track pruning
        function LMB = trackPruning(this, oldLMB)
            Len = length(oldLMB);
            r = zeros(1,Len);
            parfor idx = 1:Len
                r(idx) = oldLMB{idx}.r;
            end
            [r,index] = sort(r,'descend');
            index1 = find(r >= this.THETA);
            index = index(index1);
            Len = length(index);
            LMB = cell(1,Len);
            parfor idx = 1:Len
                LMB{idx} = oldLMB{index(idx)}
            end
        end
        % Obtain the estimated number of targets
        function N = getN(this, LMB)
            LMBLen = length(LMB);
            N = zeros(1,LMBLen);
            parfor j = 1:LMBLen
                if LMB{j}.r > 0.5
                    N(j) = 1;
                else
                    N(j) = 0;
                end
            end
            N = sum(N);
        end
    end
    
end

