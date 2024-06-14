classdef SensorCMIL < handle
    % The implementation of the fusion algorithm proposed in Reference:
    % Reference:
    % Lin Gao, Giorgio Battistelli, Luigi Chisci, "Fusion of Labeled RFS
    % Densities with Different Fields of View", IEEE Transactions on
    % aerospace and electronic systems, vol. 58, no. 6, pp. 5908-5924,
    % 2022.
    % Single sensor part
    
    properties
        MeasurementModel;
        KinematicModel;
        Pd;
        Lambda;
        ClutterModel;
        Range;
        SensorID;
        Neighbour;          % Not include itself
        BernoulliSet;
        NewbornSet;
        Ps;
        Pg;
        Gamma;
        THETA;
        N;
        CurrentGLMB;
        rBMax;
        LambdaB;
        IndexCount;
        FusionWeights;
        TD;
    end
    
    methods
        % Constructor
        function [td] = SensorCMIL(measurementModel, kinematicModel, pd, lambda, clutterModel, range, ...
                sensorID, neighbour, ps, pg, gamma, theta, n, rbmax, lambdab, fusionWeights, Td)
            td.MeasurementModel = measurementModel;
            td.KinematicModel = kinematicModel;
            td.Pd = pd;
            td.Lambda = lambda;
            td.ClutterModel = clutterModel;
            td.Range = range;
            td.SensorID = sensorID;
            td.Neighbour = neighbour;
            td.BernoulliSet = cell(0, 0);
            td.NewbornSet = cell(0, 0);
            td.Ps = ps;
            td.Pg = pg;
            td.Gamma = gamma;
            td.THETA = theta;
            td.N = n;
            td.CurrentGLMB = cell(0,0);
            td.rBMax = rbmax;
            td.LambdaB = lambdab;
            td.IndexCount = 1;
            td.FusionWeights = fusionWeights;
            td.TD = Td;
            td.setPath();
        end
        function setPath(this)
            addpath('_common\');
        end
        % Initialization
        function init(this)
            this.BernoulliSet = cell(0, 0);
            this.NewbornSet = cell(0, 0);
            this.CurrentGLMB = cell(0, 0);
            this.IndexCount = 1;
        end
        % Time update
        function [preBernoulliSet] = timeUpdate(this, Time)
            % Prediction
            oldBernoulliLen = length(this.BernoulliSet);
            NewbornLen = length(this.NewbornSet);
            preBernoulliSetOld = cell(1, oldBernoulliLen);
            preBernoulliSetNew = cell(1, NewbornLen);
            % prediction for existed targets
            parfor idx = 1:oldBernoulliLen
                preBernoulliSetOld{idx} = this.updateBernoulliComponent(this.BernoulliSet{idx});
            end
            % prediction for newborn targets
            parfor idx = 1:NewbornLen
                preBernoulliSetNew{idx} = this.genNewbornTarget(idx, Time);
            end
            this.IndexCount = NewbornLen + 1;
            preBernoulliSet = this.combineBernoulliSet(preBernoulliSetOld, preBernoulliSetNew);
        end
        % Measurement update
        function [GLMB, newbornSet] = measurementUpdate(this, preBernoulliSet, measure)
            preGLMB = this.LMB2GLMB(preBernoulliSet, this.N);
            [GLMB, newbornSet] = this.getUpdateGLMB(preGLMB, measure);
            this.CurrentGLMB = GLMB;
            this.NewbornSet = newbornSet;
        end
        % Fusion method
        function [fusedGLMB] = fusion(this, GLMBset, AdaptiveWeightFlag, Time)
            AgentNum = length(GLMBset);
            BaseDistribution = this.CurrentGLMB;
            w1 = this.FusionWeights(1);
            for idx = 1:AgentNum
                w2 = this.FusionWeights(idx + 1);
                BaseDistribution = this.pairwiseFusion(BaseDistribution, GLMBset{idx}, w1, w2, idx, AdaptiveWeightFlag, Time);
                w1 = w1 + w2;
            end
            fusedGLMB = BaseDistribution;
%             this.CurrentGLMB = fusedGLMB;
        end
        function [LMB] = GLMB2LMB(this, GLMB)
            % Normalization coefficient
            wSum = 0;
            for i = 1:length(GLMB)
                tempGLMB = GLMB{i};
                wSum = wSum + tempGLMB.w;
            end
            % Grouping
            LMBgroup = cell(0,0);
            for i = 1:length(GLMB)
                tempGLMB = GLMB{i};
                refLabelSet = tempGLMB.Label;
                refLabelSetLen = size(refLabelSet, 1);
                for n = 1:refLabelSetLen
                    LMBgroupLen = length(LMBgroup);
                    flag = false;
                    for p = 1:LMBgroupLen
                        if isequal(refLabelSet(n,:), LMBgroup{p}.l)
                            BernoulliComponentLen = length(LMBgroup{p}.BernoulliComponent);
                            LMBgroup{p}.BernoulliComponent{BernoulliComponentLen + 1}.x = tempGLMB.Distribution{n}.x;
                            LMBgroup{p}.BernoulliComponent{BernoulliComponentLen + 1}.P = tempGLMB.Distribution{n}.P;
                            LMBgroup{p}.BernoulliComponent{BernoulliComponentLen + 1}.r = tempGLMB.w / wSum;
%                             LMBgroup{p}.BernoulliComponent{BernoulliComponentLen + 1}.r = ...
%                                 double(vpa(LMBgroup{p}.BernoulliComponent{BernoulliComponentLen + 1}.r, 7));
                            flag = true;
                            break;
                        end
                    end
                    if ~flag
                        LMBgroup{LMBgroupLen + 1}.l = refLabelSet(n,:);
                        LMBgroup{LMBgroupLen + 1}.BernoulliComponent = cell(1, 1);
                        LMBgroup{LMBgroupLen + 1}.BernoulliComponent{1}.x = tempGLMB.Distribution{n}.x;
                        LMBgroup{LMBgroupLen + 1}.BernoulliComponent{1}.P = tempGLMB.Distribution{n}.P;
                        LMBgroup{LMBgroupLen + 1}.BernoulliComponent{1}.r = tempGLMB.w / wSum;
%                         LMBgroup{LMBgroupLen + 1}.BernoulliComponent{1}.r = ...
%                             double(vpa(LMBgroup{LMBgroupLen + 1}.BernoulliComponent{1}.r, 7));
                    end
                end
            end
            % Merging
            LMBgroupLen = length(LMBgroup);
            LMB = cell(1, LMBgroupLen);
            parfor idx = 1:LMBgroupLen
                x = 0;
                P = 0;
                r = 0;
                bernoulliSet = LMBgroup{idx}.BernoulliComponent;
                bernoulliSetLen = length(bernoulliSet);
                for i = 1:bernoulliSetLen
                    x = x + bernoulliSet{i}.r * bernoulliSet{i}.x;
                    r = r + bernoulliSet{i}.r;
                end
                r = real(r);
                if r > 0
                    x = x / r;
                    for i = 1:bernoulliSetLen
                        P = P + bernoulliSet{i}.r * (bernoulliSet{i}.P + (bernoulliSet{i}.x - x) * (bernoulliSet{i}.x - x)');
                    end
                    P = P / r;
                    LMB{idx}.x = x;
                    LMB{idx}.P = P;
                    LMB{idx}.l = LMBgroup{idx}.l;
                    if r > 1
                        r = 1-eps;
                    end
                    LMB{idx}.r = r;
                    LMB{idx}.r = double(vpa(LMB{idx}.r, 7));
                else
                    LMB{idx}.x = [];
                    LMB{idx}.P = [];
                    LMB{idx}.l = LMBgroup{idx}.l;
                    LMB{idx}.r = 0;
                end
            end
            % 剪枝
%             LMB = this.trackPruning(LMB);
            this.BernoulliSet = LMB;
        end
        function [transportDataGLMB] = genTransportDataGLMB(this)
            transportDataGLMB = this.CurrentGLMB;
        end
        function neighbour = getNeighbour(this)
            neighbour = this.Neighbour;
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
                LMB{idx} = oldLMB{index(idx)};
                if LMB{idx}.r >= 1
                    LMB{idx}.r = 1-eps;
                end
%                 LMB{idx}.r = vpa(LMB{idx}.r, 7);
            end
        end
    end
    
    methods
        %% Time Update
        % Prediction for existed targets
        function [bernoulliComponent] = updateBernoulliComponent(this, oldBernoulliComponent)
            bernoulliComponent.r = this.Ps * oldBernoulliComponent.r;
            bernoulliComponent.l = oldBernoulliComponent.l;
            bernoulliComponent.x = this.KinematicModel.estimateNewState(oldBernoulliComponent.x);
            F = this.KinematicModel.getLinear(oldBernoulliComponent.x);
            bernoulliComponent.P = F * oldBernoulliComponent.P * F' + this.KinematicModel.CovMatrix;
        end
        % Prediction for newborn targets
        function [bernoulliComponent] = genNewbornTarget(this, index, time)
            bernoulliComponent.x = this.NewbornSet{index}.x;
            bernoulliComponent.P = this.NewbornSet{index}.P;
            bernoulliComponent.r = this.NewbornSet{index}.r;
            bernoulliComponent.l = [time, index];
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
        %% Measurement Update
        % LMB RFS to GLMB RFS
        function [GLMB] = LMB2GLMB(this, BernoulliSet, Num)
            BernoulliSetLen = length(BernoulliSet);
            r = zeros(BernoulliSetLen,1);
            Distribution = cell(1,BernoulliSetLen);
            Label = zeros(BernoulliSetLen,2);
            for idx = 1:BernoulliSetLen
                Distribution{idx}.x = BernoulliSet{idx}.x;
                Distribution{idx}.P = BernoulliSet{idx}.P;
                Label(idx,:) = BernoulliSet{idx}.l;
                r(idx) = BernoulliSet{idx}.r;
            end
            GLMB.Distribution = Distribution;
            GLMB.l = Label;
            cost = r./(1-r);
            neglogcost = -log(cost);
            [paths, nlcost] = kshortestwrap_pred(neglogcost, Num);
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
        % Get the posterior GLMB distribution
        function [GLMB, newbornSet] = getUpdateGLMB(this, preGLMB, measure)
            measureLen = size(measure,2);
            % update each Bernoulli Component
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
                Z = [Z, ones(1,BernoulliLen) * idx];
            end
            % Update the GLMB weight
            I_len = length(preGLMB.I);
            posteriorGLMB = cell(1,I_len);
            parfor idx = 1:I_len
                posteriorGLMB{idx} = this.getPosteriorGLMB(preGLMB.l(preGLMB.I{idx},:), preGLMB.w(idx), Distribution, Label, ETA, Z, 5);
            end
            % 归一化
            wSum = 0;
            for idx = 1:I_len
                wSum = wSum + sum(posteriorGLMB{idx}.w);
            end
            parfor idx = 1:I_len
                posteriorGLMB{idx}.w = posteriorGLMB{idx}.w / wSum;
                posteriorGLMB{idx}.w = double(vpa(posteriorGLMB{idx}.w, 7));
            end
            newbornSet = this.genBirthDistribution(posteriorGLMB, measure);
            GLMB = this.simplifiedGLMB(posteriorGLMB, Distribution, Label, Z);
            updateGLMB = this.pruningGLMB(GLMB);
            GLMB = updateGLMB;
        end
        % Single target update 2
        function [Distribution, Label, ETA] = updateBernoulli2(this, GLMB, measure, index, BernoulliLen)
            Distribution = cell(1,BernoulliLen);
            Label = zeros(BernoulliLen,2);
            ETA = zeros(1, BernoulliLen);
            for idx2 = 1:BernoulliLen
                [Distribution{idx2}, eta] = this.updateBernoulli1(GLMB.Distribution{idx2}, measure, index);
                Label(idx2,:) = GLMB.l(idx2,:);
                ETA(idx2) = eta;
            end
        end
        % Single target update 1
        function [newDistribution, eta] = updateBernoulli1(this, oldDistribution, measure, measureIdx)
            if measureIdx == 0
                eta = 1 - this.getPd(oldDistribution.x);
                newDistribution = oldDistribution;
            else
                P_pre = oldDistribution.P;
                x_pre = oldDistribution.x;
                H = this.MeasurementModel.getLinear(x_pre);
                R = this.MeasurementModel.CovMatrix;
                S = H*P_pre*H' + R;
                S = (S + S') / 2;
                Vs = chol(S);
                inv_sqrt_S = inv(Vs);
                iS= inv_sqrt_S*inv_sqrt_S';
                K = P_pre * H' * iS;
                z = this.MeasurementModel.estimateNewState(x_pre);
                dz = measure(:,measureIdx) - z;
                x_est = x_pre + K * dz;
                P_est = P_pre - K * H * P_pre;
                newDistribution.x = x_est;
                newDistribution.P = P_est;
                temp = 1/sqrt(det(2*pi*S)) * exp(-1/2 * dz' * iS * dz) + eps;
                eta = this.getPd(x_pre) * temp / (this.ClutterModel.getProb(z) * this.Lambda);
            end
        end
        function [pd] = getPd(this, x)
            z = this.MeasurementModel.estimateNewState(x);
            zLen = length(z);
            flag = true;
            for idx = 1:zLen
                if z(idx) < this.Range(idx,1) || z(idx) > this.Range(idx,2)
                    flag = false;
                    break;
                end
            end
            if flag
                pd = this.Pd;
            else
                pd = 0;
            end
        end
        % Find the GLMB components with the Num largest weights.
        function [GLMB] = getPosteriorGLMB(this, LabelSet, w, Distribution, Label, ETA, Z, Num)
            % LabelSet: The label set for a specific GLMB component
            GLMB.Distribution = Distribution;
            GLMB.Label = Label;
            LabelSetLen = size(LabelSet, 1);
            if LabelSetLen == 0
                % No target assumption
                GLMB.w = w;
                GLMB.I = {[]};
                GLMB.Z = [];
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
                [assignments, cost] = murty_custom(C, Num);
                costLen = length(cost);
                newPath = cell(1,costLen);
                tgtNum = size(assignments,2);
                parfor j = 1:costLen
                    tmpIndexList = assignments(j,:);
                    P = [];
                    for k = 1:tgtNum
                        if tmpIndexList(k) > measureNum
                            node = [k;0];
                        else
                            node = [k; tmpIndexList(k)];
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
                GLMB.Z = Z;
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
        % Obtain the newborn density
        function [newbornSet] = genBirthDistribution(this, GLMB, measure)
            measureLen = size(measure,2);
            rU = zeros(1,measureLen);
            GLMBLen = length(GLMB);
            for i = 1:measureLen
                temp = 0;
                for j = 1:GLMBLen
                    ILen = length(GLMB{j});
                    for k = 1:ILen
                        if ~isempty(find(GLMB{j}.I{k} == i))
                            temp = temp + GLMB{j}.w(k);
                        end
                    end
                end
                rU(i) = temp;
            end
            rB = 1-rU;
            rB = rB / sum(rB);
            newbornSet = cell(1, measureLen);
            for idx = 1:measureLen
                newbornSet{idx}.x = this.genSingleTargetNewbornState(measure(:,idx));
                newbornSet{idx}.P = diag([100; 100; 100; 100].^2);
                newbornSet{idx}.r = min(this.rBMax, rB(idx) * this.LambdaB);
                newbornSet{idx}.r = double(vpa(newbornSet{idx}.r, 7));
            end
            % 剪枝
            tempNewbornSet = cell(0,0);
            count = 1;
            for idx = 1:measureLen
                if newbornSet{idx}.r > 0.1
                    tempNewbornSet{count} = newbornSet{idx};
                    count = count + 1;
                end
            end
            newbornSet = tempNewbornSet;
            this.NewbornSet = newbornSet;
        end
        function x = genSingleTargetNewbornState(this, z)
            pos = this.MeasurementModel.getCartesianPos(z);
            x = [];
            dimension = length(pos);
            for idx = 1:dimension
                tempX = [pos(idx);0];
                x = [x; tempX];
            end
        end
        % Simplified the GLMB components
        function GLMB = simplifiedGLMB(this, oldGLMB, Distribution, Label, Z)
            Len = length(oldGLMB);
            GLMB = cell(1,Len);
            wSum = 0;
            for idx = 1:Len
                GLMBcomponent = oldGLMB{idx};
                wLen = length(GLMBcomponent.w);
                LMBgroup = cell(0,0);
                w = 0;
                % Grouping
                for i = 1:wLen
                    tempLabelSet = Label(GLMBcomponent.I{i},:);
                    tempLabelSetLen = size(tempLabelSet,1);
                    LMBgroupLen = length(LMBgroup);
                    w = w + GLMBcomponent.w(i);
                    tempI = GLMBcomponent.I{i};
                    for j = 1:tempLabelSetLen
                        flag = true;
                        for k = 1:LMBgroupLen
                            if isequal(tempLabelSet(j,:), LMBgroup{k}.Label)
                                flag = false;
                                tempLen = length(LMBgroup{k}.BernoulliComponent);
                                LMBgroup{k}.BernoulliComponent{tempLen + 1}.x = GLMBcomponent.Distribution{tempI(j)}.x;
                                LMBgroup{k}.BernoulliComponent{tempLen + 1}.P = GLMBcomponent.Distribution{tempI(j)}.P;
                                LMBgroup{k}.BernoulliComponent{tempLen + 1}.r = GLMBcomponent.w(i);
                                LMBgroup{k}.BernoulliComponent{tempLen + 1}.z = GLMBcomponent.Z(tempI(j));
                                break;
                            end
                        end
                        if flag
                            tempLMBgroupLen = length(LMBgroup);
                            LMBgroup{tempLMBgroupLen + 1}.Label = tempLabelSet(j,:);
                            LMBgroup{tempLMBgroupLen + 1}.BernoulliComponent = cell(0,0);
                            LMBgroup{tempLMBgroupLen + 1}.BernoulliComponent{1}.x = GLMBcomponent.Distribution{tempI(j)}.x;
                            LMBgroup{tempLMBgroupLen + 1}.BernoulliComponent{1}.P = GLMBcomponent.Distribution{tempI(j)}.P;
                            LMBgroup{tempLMBgroupLen + 1}.BernoulliComponent{1}.r = GLMBcomponent.w(i);
                            LMBgroup{tempLMBgroupLen + 1}.BernoulliComponent{1}.z = GLMBcomponent.Z(tempI(j));
                        end
                    end
                end
                wSum = wSum + w;
                % Merging
                GLMB{idx}.w = w;
                if w > 0
                    newLabelSet = [];
                    newDistribution = cell(0,0);
                    newZ = cell(0,0);
                    LMBgroupLen = length(LMBgroup);
                    for i = 1:LMBgroupLen
                        x = 0;
                        P = 0;
                        Z = [];
                        newLabelSet = [newLabelSet; LMBgroup{i}.Label];
                        BernoulliComponent = LMBgroup{i}.BernoulliComponent;
                        BernoulliComponentLen = length(BernoulliComponent);
                        for j = 1:BernoulliComponentLen
                            x = x + BernoulliComponent{j}.r * BernoulliComponent{j}.x;
                            Z = [Z, BernoulliComponent{j}.z];
                        end
                        x = x / w;
                        for j = 1:BernoulliComponentLen
                            P = P + BernoulliComponent{j}.r * ( BernoulliComponent{j}.P + (BernoulliComponent{j}.x - x) * (BernoulliComponent{j}.x - x)' );
                        end
                        P = P / w;
                        newDistribution{i}.x = x;
                        newDistribution{i}.P = P;
                        newZ{i} = Z;
                    end
                    GLMB{idx}.Label = newLabelSet;
                    GLMB{idx}.Distribution = newDistribution;
                    GLMB{idx}.Z = newZ;
                else
                    tempI = GLMBcomponent.I{1};
                    GLMB{idx}.Label = Label(tempI,:);
                    GLMB{idx}.Distribution = cell(1,length(tempI));
                    GLMB{idx}.Z = cell(1, length(tempI));
                    for k = 1:length(tempI)
                        GLMB{idx}.Distribution{k} = Distribution{tempI(k)};
                        try
                            GLMB{idx}.Z{k} = Z(tempI(k));
                        catch
                            disp('Stop!');
                        end
                    end
                end
            end
            GLMBLen = length(GLMB);
            for idx = 1:GLMBLen
                GLMB{idx}.w = GLMB{idx}.w / wSum;
                GLMB{idx}.w = double(vpa(GLMB{idx}.w, 7));
            end
        end
        function newGLMB = pruningGLMB(this, GLMB)
            GLMBLen = length(GLMB);
            newGLMB = cell(0,0);
            wSum = 0;
            count = 1;
            for idx = 1:GLMBLen
                if GLMB{idx}.w > 0.001
                    newGLMB{count} = GLMB{idx};
                    wSum = wSum + GLMB{idx}.w;
                    count = count + 1;
                end
            end
            newGLMBLen = length(newGLMB);
            for idx = 1:newGLMBLen
                newGLMB{idx}.w = newGLMB{idx}.w / wSum;
            end
        end
        %% Fusion
        function [fusedGLMB] = pairwiseFusion(this, BaseGLMB, AppendGLMB, w1, w2, agent, AdaptiveWeightFlag, Time)
            % Match labels between two LRFS densities by solving the
            % assignment problem.
            wSum = w1 + w2;
            w1 = w1 / wSum;
            w2 = w2 / wSum;
            BaseLMB = this.GLMB2LMB(BaseGLMB);
            AppendLMB = this.GLMB2LMB(AppendGLMB);
            BaseLMBLen = length(BaseLMB);
            AppendLMBLen = length(AppendLMB);
            [AssignmentMatrix, C_star] = this.optimalLabelMatching(BaseLMB, AppendLMB);
            %%%%%%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             figure(20 + agent);
%             cla;
%             hold on
%             grid on
%             for idx = 1:BaseLMBLen
%                 plot(BaseLMB{idx}.x(1), BaseLMB{idx}.x(3),'Marker','*', 'Color','r');
%                 txt = sprintf('(%d,%d),r=%.2f', BaseLMB{idx}.l(1), BaseLMB{idx}.l(2),BaseLMB{idx}.r);
%                 text(BaseLMB{idx}.x(1), BaseLMB{idx}.x(3), txt);
%             end
%             for idx = 1:AppendLMBLen
%                 plot(AppendLMB{idx}.x(1), AppendLMB{idx}.x(3), 'Marker', 'o', 'Color', 'b');
%                 txt = sprintf('(%d,%d),r=%.2f', AppendLMB{idx}.l(1), AppendLMB{idx}.l(2), AppendLMB{idx}.r);
%                 text(AppendLMB{idx}.x(1), AppendLMB{idx}.x(3), txt);
%             end
%             hold off
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Find the disjoint label subspaces L_m (for m=1,2,3) according
            % to the results of label matching
            LabelSpaceE1 = cell(0,0);
            LabelSpaceE2 = cell(0,0);
            LabelSpaceC = cell(0,2);
            % Construct the lookup tables
            TableE1 = [];
            TableE2 = [];
            TableC1 = [];
            TableC2 = [];
            for idx1 = 1:(BaseLMBLen + 1)
                for idx2 = 1:(AppendLMBLen + 1)
                    if AssignmentMatrix(idx1, idx2) == 1
                        if idx1 == BaseLMBLen + 1
                            % Labels in the exclusive label space of
                            % neighbour
                            LabelSpaceE2Len = length(LabelSpaceE2);
                            LabelSpaceE2{LabelSpaceE2Len + 1} = AppendLMB{idx2}.l;
                            tmpL = AppendLMB{idx2}.l;
                            TableE2(tmpL(1), tmpL(2)) = LabelSpaceE2Len + 1;
                            TableC2(tmpL(1), tmpL(2)) = 0;
                        elseif idx2 == AppendLMBLen + 1
                            % Labels in the exclusive label space of the
                            % local agent
                            LabelSpaceE1Len = length(LabelSpaceE1);
                            LabelSpaceE1{LabelSpaceE1Len + 1} = BaseLMB{idx1}.l;
                            tmpL = BaseLMB{idx1}.l;
                            TableE1(tmpL(1), tmpL(2)) = LabelSpaceE1Len + 1;
                            TableC1(tmpL(1), tmpL(2)) = 0;
                        else
                            % Assigned labels
                            [Num, SpaceNum] = size(LabelSpaceC);
                            LabelSpaceC{Num+1, 1} = BaseLMB{idx1}.l;
                            LabelSpaceC{Num+1, 2} = AppendLMB{idx2}.l;
                            tmpL1 = BaseLMB{idx1}.l;
                            tmpL2 = AppendLMB{idx2}.l;
                            TableE1(tmpL1(1), tmpL1(2)) = 0;
                            TableC1(tmpL1(1), tmpL1(2)) = Num+1;
                            TableE2(tmpL2(1), tmpL2(2)) = 0;
                            TableC2(tmpL2(1), tmpL2(2)) = Num+1;
                        end
                    end
                end
            end
            % Decompose \pi^i and \pi^j to \pi^i={\pi_m^i}_{m=1}^3 and
            % \pi^j={\pi_m^j}_{m=1}^3
            [BaseGLMBE, BaseGLMBC] = this.decomposeGLMB(BaseGLMB, TableE1, TableC1);
            [AppendGLMBE, AppendGLMBC] = this.decomposeGLMB(AppendGLMB, TableE2, TableC2);
            % Fuse each local subdensity pair (\pi_m^i, \pi_m^j), for
            % m=1,2,3, into the global subdensity \bar{\pi_m}
            fusedGLMB1 = cell(0,0);
            count = 1;
            for idx = 1:length(BaseGLMBE)
                if ~isempty(BaseGLMBE{idx}.Label)
                    if AdaptiveWeightFlag
                        fusedGLMB1{count} = BaseGLMBE{idx};
                        count = count + 1;
                    else
                        fusedGLMB1{count} = BaseGLMBE{idx};
                        fusedGLMB1{count}.w = fusedGLMB1{count} * w1;
                        count = count + 1;
                    end
                end
            end
            w_empty = 1;
            fusedGLMB1Len = length(fusedGLMB1);
            for idx = 1:fusedGLMB1Len
                w_empty = w_empty - fusedGLMB1{idx}.w;
            end
            if w_empty < 0
                w_empty = 0;
                fusedGLMB1{fusedGLMB1Len + 1}.w = w_empty;
                fusedGLMB1{fusedGLMB1Len + 1}.Label = [];
                fusedGLMB1{fusedGLMB1Len + 1}.Distribution = cell(0,0);
                % Normalization 
                wSum = 0;
                for idx = 1:fusedGLMB1Len
                    wSum = wSum + fusedGLMB1{idx}.w;
                end
                for idx = 1:fusedGLMB1Len
                    fusedGLMB1{idx}.w = fusedGLMB1{idx}.w / wSum; 
                end
            else
                fusedGLMB1{fusedGLMB1Len + 1}.w = w_empty;
                fusedGLMB1{fusedGLMB1Len + 1}.Label = [];
                fusedGLMB1{fusedGLMB1Len + 1}.Distribution = cell(0,0);
            end
            fusedGLMB1 = this.GLMBdataValidation(fusedGLMB1);
            fusedGLMB2 = cell(0,0);
            count = 1;
            AppendLabelSpace = cell(1, length(LabelSpaceE2));
            TableAppend = [];
            % Generating AppendLabelSet and TableAppend
            for idx = 1:length(LabelSpaceE2)
                % IndexCount
                tmpLabel = [Time, this.IndexCount];
                AppendLabelSpace{idx} = tmpLabel;
                TableAppend(tmpLabel(1), tmpLabel(2)) = idx;
                this.IndexCount = this.IndexCount + 1;
            end
            for idx = 1:length(AppendGLMBE)
                if ~isempty(AppendGLMBE{idx}.Label)
                    tmpLabel = [];
                    Number = size(AppendGLMBE{idx}.Label, 1);
                    for i = 1:Number
                        tmpIndex = TableE2(AppendGLMBE{idx}.Label(i,1), AppendGLMBE{idx}.Label(i,2));
                        if tmpIndex == 0
                            disp('There are errors when generating TableE2!');
                        end
                        tmpLabel = [tmpLabel; AppendLabelSpace{tmpIndex}];
                    end
                    fusedGLMB2{count}.Label = tmpLabel;
                    fusedGLMB2{count}.Distribution = AppendGLMBE{idx}.Distribution;
                    if AdaptiveWeightFlag
                        fusedGLMB2{count}.w = AppendGLMBE{idx}.w;
                    else
                        fusedGLMB2{count}.w = w2 * AppendGLMBE{idx}.w;
                    end
                    count = count + 1;
                end
            end
            fusedGLMB2Len = length(fusedGLMB2);
            w_empty = 1;
            for idx = 1:fusedGLMB2Len
                w_empty = w_empty - fusedGLMB2{idx}.w;
            end
            if w_empty < 0
                w_empty = 0;
                fusedGLMB2{fusedGLMB2Len + 1}.w = w_empty;
                fusedGLMB2{fusedGLMB2Len + 1}.Label = [];
                fusedGLMB2{fusedGLMB2Len + 1}.Distribution = cell(0,0);
                % Normalization 
                wSum = 0;
                for idx = 1:fusedGLMB2Len
                    wSum = wSum + fusedGLMB2{idx}.w;
                end
                for idx = 1:fusedGLMB2Len
                    fusedGLMB2{idx}.w = fusedGLMB2{idx}.w / wSum; 
                end
            else
                fusedGLMB2{fusedGLMB2Len + 1}.w = w_empty;
                fusedGLMB2{fusedGLMB2Len + 1}.Label = [];
                fusedGLMB2{fusedGLMB2Len + 1}.Distribution = cell(0,0);
            end
            fusedGLMB2 = this.GLMBdataValidation(fusedGLMB2);
            fusedCommonGLMB = this.fusionCommonGLMB(BaseGLMBC, AppendGLMBC, TableC1, TableC2, w1, w2, LabelSpaceC);
%             fusedCommonGLMBLen = length(fusedCommonGLMB);
%             FusedGLMBLen = length(fusedGLMB);
%             for idx = 1:fusedCommonGLMBLen
%                 fusedGLMB{FusedGLMBLen + idx} = fusedCommonGLMB{idx};
%             end
            % Obtain the fused density as \bar{\pi}=\prod_{m=1}^3\bar{\pi_m}
%             this.checkData(fusedGLMB1);
%             this.checkData(fusedGLMB2);
%             this.checkData(fusedCommonGLMB);
            [fusedGLMB] = this.multiplyGLMB(fusedGLMB1, fusedGLMB2);
            [fusedGLMB] = this.multiplyGLMB(fusedGLMB, fusedCommonGLMB);
            % Normalization
            fusedGLMB = this.GLMBdataValidation(fusedGLMB);
            
        end
        % Solve the label assignment problem
        function [AssignmentMatrix, CostMatrix] = optimalLabelMatching(this, baseLMB, appendLMB)
            baseLMBLen = length(baseLMB);
            appendLMBLen = length(appendLMB);
            CostMatrix = zeros(baseLMBLen + 1, appendLMBLen + 1);
            for idx1 = 1:(baseLMBLen + 1)
                for idx2 = 1:(appendLMBLen + 1)
                    if idx1 <= baseLMBLen && idx2 <= appendLMBLen
                        x1 = baseLMB{idx1}.x;
                        P1 = baseLMB{idx1}.P;
                        x2 = appendLMB{idx2}.x;
                        P2 = appendLMB{idx2}.P;
                        % Cauchy-Schwarz divergence
                        zmk = 1/sqrt(det(2*pi*(P1 + P2))) * exp(-1/2 * (x1-x2)' * pinv(P1 + P2) * (x1-x2));
                        CSD = -log(zmk) + 1/2 * log(sqrt(det(pinv(P1) / (2*pi)))) + 1/2 * log(sqrt(det(pinv(P2) / (2*pi))));
                        CostMatrix(idx1, idx2) = CSD;
                    else
                        if idx1 > baseLMBLen && idx2 <= appendLMBLen
                            CostMatrix(idx1, idx2) = this.TD;
                        elseif idx1 <= baseLMBLen && idx2 > appendLMBLen
                            CostMatrix(idx1, idx2) = this.TD;
                        else
                            CostMatrix(idx1, idx2) = inf;
                        end
                    end
                end
            end
            % Obtain the label assignment according to the Hungarian
            % algorithm.
            % Convert the CostMatrix to another matrix that can be
            % processed by the Hungarian algorithm directly.
            A = CostMatrix(1:baseLMBLen, 1:appendLMBLen);
            B1 = CostMatrix(1:baseLMBLen, appendLMBLen + 1);
            B1 = -log(eye(baseLMBLen)) + diag(B1);
            B2 = CostMatrix(baseLMBLen + 1, 1:appendLMBLen);
            B2 = -log(eye(appendLMBLen)) + diag(B2);
            C = ones(appendLMBLen, baseLMBLen) * this.TD;
            newCostMatrix = [A, B1; B2, C];
            [~, x] = HungarianAlgorithm(newCostMatrix);
            xA = x(1:baseLMBLen, 1:appendLMBLen);
            xB1 = x(1:baseLMBLen, (appendLMBLen+1):end);
            xB1 = sum(xB1,2);
            xB2 = x((baseLMBLen + 1):end, 1:appendLMBLen);
            xB2 = sum(xB2,1);
            AssignmentMatrix = [xA, xB1; xB2, 0];
        end
        % Decompose the GLMB RFS
        function [GLMBE, GLMBC] = decomposeGLMB(this, GLMB, TableE, TableC)
            GLMBEset = cell(0,0);
            GLMBCset = cell(0,0);
            GLMBLen = length(GLMB);
            % GLMB component: w, Label, Distribution, Z(neglected)
            for idx = 1:GLMBLen
                tmpGLMB = GLMB{idx};
                tmpLabel = tmpGLMB.Label;
                tmpLabelE = [];
                tmpLabelC = [];
                tmpIndexE = [];
                tmpIndexC = [];
                DistributionE = cell(0,0);
                DistributionC = cell(0,0);
                [Len,~] = size(tmpLabel);
                for i = 1:Len
                    l = tmpLabel(i,:);
                    refIdx1 = TableE(l(1), l(2));
                    refIdx2 = TableC(l(1), l(2));
                    if refIdx2 > 0
                        % The label belongs to the common label space
                        tmpLabelC = [tmpLabelC; l];
                        tmpIndexC = [tmpIndexC, refIdx2];
                        DistributionCLen = length(DistributionC);
                        DistributionC{DistributionCLen + 1} = tmpGLMB.Distribution{i};
                    elseif refIdx1 > 0
                        % The label belongs to the exclusive label space
                        tmpLabelE = [tmpLabelE; l];
                        tmpIndexE = [tmpIndexE, refIdx1];
                        DistributionELen = length(DistributionE);
                        DistributionE{DistributionELen + 1} = tmpGLMB.Distribution{i};
                    else
                        disp('Not likely to happen!');
                    end
                end
                % Append to the sets
                [GLMBEset] = this.append2GLMBsets(tmpGLMB.w, GLMBEset, tmpLabelE, tmpIndexE, DistributionE);
                [GLMBCset] = this.append2GLMBsets(tmpGLMB.w, GLMBCset, tmpLabelC, tmpIndexC, DistributionC);
            end
            % Reorganizing
            GLMBE = this.reorganizingDecomposedGLMB(GLMBEset);
            GLMBC = this.reorganizingDecomposedGLMB(GLMBCset);
        end
        function [GLMBsets] = append2GLMBsets(this, w, oldGLMBsets, Label, Index, Distribution)
            GLMBsets = oldGLMBsets;
            if ~isempty(Index)
                oldGLMBsetsLen = length(oldGLMBsets);
                flag = true;
                for i = 1:oldGLMBsetsLen
                    if length(intersect(Index, oldGLMBsets{i}.IndexSet)) == length(union(Index, oldGLMBsets{i}.IndexSet))
                        % Matched
                        tmpLen = length(oldGLMBsets{i}.IndexSet);
                        ComponentsLen = length(oldGLMBsets{i}.Components);
                        GLMBsets{i}.Components{ComponentsLen + 1}.w = w;
                        GLMBsets{i}.Components{ComponentsLen + 1}.Distribution = cell(1, tmpLen);
                        for j = 1:tmpLen
                            tmpID = find(Index == GLMBsets{i}.IndexSet(j));
                            GLMBsets{i}.Components{ComponentsLen + 1}.Distribution{j} = Distribution{tmpID(1)};
                        end
                        flag = false;
                        break;
                    end
                end
                if flag
                    GLMBsets{oldGLMBsetsLen + 1}.IndexSet = Index;
                    GLMBsets{oldGLMBsetsLen + 1}.Label = Label;
                    GLMBsets{oldGLMBsetsLen + 1}.Components = cell(1,1);
                    GLMBsets{oldGLMBsetsLen + 1}.Components{1}.w = w;
                    GLMBsets{oldGLMBsetsLen + 1}.Components{1}.Distribution = Distribution;
                end
            end
        end
        function [GLMB] = reorganizingDecomposedGLMB(this, GLMBset)
            % GLMBset: IndexSet, Label, Components
            % Components: w, Distribution
            GLMBsetLen = length(GLMBset);
            GLMB = cell(1, GLMBsetLen);
            w_empty = 1;
            for idx = 1:GLMBsetLen
                ComponentsLen = length(GLMBset{idx}.Components);
                w = zeros(1, ComponentsLen);
                for i = 1:ComponentsLen
                    w(i) = GLMBset{idx}.Components{i}.w;
                end
                sumW = sum(w);
                w = w ./ sumW;
                w_empty = w_empty - sumW;
                GLMB{idx}.w = sumW;
                GLMB{idx}.Label = GLMBset{idx}.Label;
                Number = length(GLMBset{idx}.IndexSet);
                Distribution = cell(1, Number);
                parfor i = 1:Number
                    x = 0;
                    P = 0;
                    for j = 1:ComponentsLen
                        x = x + w(j) * GLMBset{idx}.Components{j}.Distribution{i}.x;
                    end
                    for j = 1:ComponentsLen
                        tmp_P = GLMBset{idx}.Components{j}.Distribution{i}.P;
                        tmp_x = GLMBset{idx}.Components{j}.Distribution{i}.x;
                        P = P + w(j) * (tmp_P + (x - tmp_x) * (x - tmp_x)');
                    end
                    Distribution{i}.x = x;
                    Distribution{i}.P = P;
                end
                GLMB{idx}.Distribution = Distribution;
            end
            GLMB{GLMBsetLen + 1}.w = w_empty;
            GLMB{GLMBsetLen + 1}.Distribution = cell(0,0);
            GLMB{GLMBsetLen + 1}.Label = [];
            GLMB = this.GLMBdataValidation(GLMB);
        end
        function [fusedGLMB] = fusionCommonGLMB(this, GLMB1, GLMB2, Table1, Table2, w1, w2, LabelSpaceC)
            GLMB1Len = length(GLMB1);
            GLMB2Len = length(GLMB2);
            % Appending IndexSet
            for idx = 1:GLMB1Len
                tmpIndexSet = [];
                Number = size(GLMB1{idx}.Label,1);
                for i = 1:Number
                    ID = Table1(GLMB1{idx}.Label(i,1), GLMB1{idx}.Label(i,2));
                    if ID == 0
                        disp('There are errors when generating Table1!');
                    end
                    tmpIndexSet = [tmpIndexSet, ID];
                end
                GLMB1{idx}.IndexSet = tmpIndexSet;
            end
            for idx = 1:GLMB2Len
                tmpIndexSet = [];
                Number = size(GLMB2{idx}.Label,1);
                for i = 1:Number
                    ID = Table2(GLMB2{idx}.Label(i,1), GLMB2{idx}.Label(i,2));
                    if ID == 0
                        disp('There are errors when generating Table2!');
                    end
                    tmpIndexSet = [tmpIndexSet, ID];
                end
                GLMB2{idx}.IndexSet = tmpIndexSet;
            end
            % GLMB: w, Distribution, Label, IndexSet
            usedFlag = zeros(1, GLMB2Len);
            fusedGLMB = cell(0,0);
            count = 1;
            for idx1 = 1:GLMB1Len
                tmpGLMB1 = GLMB1{idx1};
                Flag = 1;
                for idx2 = 1:GLMB2Len
                    if usedFlag(idx2) == 0
                        % unmatched GLMB components in GLMB2
                        tmpGLMB2 = GLMB2{idx2};
                        if length(intersect(tmpGLMB1.IndexSet, tmpGLMB2.IndexSet)) == length(union(tmpGLMB1.IndexSet, tmpGLMB2.IndexSet)) &&...
                                ~isempty(tmpGLMB1.IndexSet)
                            % matched
                            wFused = tmpGLMB1.w * w1 + tmpGLMB2.w * w2;
                            weights = [tmpGLMB1.w * w1 / wFused, tmpGLMB2.w * w2 / wFused];
                            Number = length(tmpGLMB1.IndexSet);
                            Distribution = cell(1, Number);
                            for i = 1:Number
                                x1 = tmpGLMB1.Distribution{i}.x;
                                P1 = tmpGLMB1.Distribution{i}.P;
                                tmpIDX = find(tmpGLMB2.IndexSet == tmpGLMB1.IndexSet(i));
                                x2 = tmpGLMB2.Distribution{tmpIDX(1)}.x;
                                P2 = tmpGLMB2.Distribution{tmpIDX(1)}.P;
                                x = weights(1) * x1 + weights(2) * x2;
                                P = weights(1) * (P1 + (x - x1) * (x - x1)') + weights(2) * (P2 + (x - x2) * (x - x2)');
                                Distribution{i}.x = x;
                                Distribution{i}.P = P;
                            end
                            fusedGLMBcomponent.w = wFused;
                            fusedGLMBcomponent.Label = tmpGLMB1.Label;
                            fusedGLMBcomponent.Distribution = Distribution;
                            fusedGLMB{count} = fusedGLMBcomponent;
                            count = count + 1;
                            Flag = 0;
                            usedFlag(idx2) = 1;
                            break;
                        end
                    end
                end
                if Flag == 1
                    fusedGLMBcomponent.w = tmpGLMB1.w * w1;
                    fusedGLMBcomponent.Label = tmpGLMB1.Label;
                    fusedGLMBcomponent.Distribution = tmpGLMB1.Distribution;
                    fusedGLMB{count} = fusedGLMBcomponent;
                    count = count + 1;
                end
            end
            count = GLMB1Len + 1;
            for idx1 = 1:GLMB2Len
                if usedFlag(idx1) == 0 && size(GLMB2{idx1}.Label,1) > 0
                    tmpGLMB2 = GLMB2{idx1};
                    fusedGLMBcomponent.w = tmpGLMB2.w * w2;
                    tmpLabel = [];
                    for i = 1:size(tmpGLMB2.Label,1)
                        tmpID = Table2(tmpGLMB2.Label(i,1), tmpGLMB2.Label(i,2));
                        tmpLabel = [tmpLabel; LabelSpaceC{tmpID,1}];
                    end
                    fusedGLMBcomponent.Label = tmpLabel;
                    fusedGLMBcomponent.Distribution = tmpGLMB2.Distribution;
                    fusedGLMB{count} = fusedGLMBcomponent;
                    count = count + 1;
                end
            end
            w_empty = 1;
            fusedGLMBLen = length(fusedGLMB);
            for idx = 1:fusedGLMBLen
                w_empty = w_empty - fusedGLMB{idx}.w;
            end
            if w_empty < 0
                w_empty = 0;
                fusedGLMB{fusedGLMBLen + 1}.w = w_empty;
                fusedGLMB{fusedGLMBLen + 1}.Label = [];
                fusedGLMB{fusedGLMBLen + 1}.Distribution = cell(0,0);
                % Normalization
                wSum = 0;
                for idx = 1:fusedGLMBLen
                    wSum = wSum + fusedGLMB{idx}.w;
                end
                for idx = 1:fusedGLMBLen
                    fusedGLMB{idx}.w = fusedGLMB{idx}.w / wSum;
                end
            else
                fusedGLMB{fusedGLMBLen + 1}.w = w_empty;
                fusedGLMB{fusedGLMBLen + 1}.Label = [];
                fusedGLMB{fusedGLMBLen + 1}.Distribution = cell(0,0);
            end
            
        end
        function [fusedGLMB] = multiplyGLMB(this, GLMB1, GLMB2)
            GLMB1Len = length(GLMB1);
            GLMB2Len = length(GLMB2);
            w = zeros(1, GLMB1Len * GLMB2Len);
            for idx1 = 1:GLMB1Len
                for idx2 = 1:GLMB2Len
                    tmpGLMB1 = GLMB1{idx1};
                    tmpGLMB2 = GLMB2{idx2};
                    tmpW = tmpGLMB1.w * tmpGLMB2.w;
                    w((idx1 - 1) * GLMB2Len + idx2) = tmpW;
                end
            end
            % Pruning
%             [~, index] = sort(w, 'descend');
%             ComponentNum = min(this.N, GLMB1Len * GLMB2Len);
            index = find(w > 5e-4);
            ComponentNum = length(index);
            fusedGLMB = cell(1, ComponentNum);
            parfor idx = 1:ComponentNum
                ID = index(idx);
                id2 = mod(ID, GLMB2Len);
                id1 = floor(ID/GLMB2Len) + 1;
                tmpW = w(ID);
                if id2 == 0
                    id2 = GLMB2Len;
                    id1 = id1 - 1;
                end
                tmpGLMB1 = GLMB1{id1};
                tmpGLMB2 = GLMB2{id2};
                tmpLabel = [tmpGLMB1.Label; tmpGLMB2.Label];
                Number1 = size(tmpGLMB1.Label,1);
                Number2 = size(tmpGLMB2.Label,1);
                tmpDistribution = cell(1, size(tmpLabel,1));
                for i = 1:Number1
                    tmpDistribution{i} = tmpGLMB1.Distribution{i};
                end
                for i = 1:Number2
                    tmpDistribution{i + Number1} = tmpGLMB2.Distribution{i};
                end
                fusedGLMB{idx}.w = tmpW;
                fusedGLMB{idx}.Label = tmpLabel;
                fusedGLMB{idx}.Distribution = tmpDistribution;
            end
        end
        function [newGLMB] = GLMBdataValidation(this, GLMB)
            GLMBLen = length(GLMB);
            newGLMB = cell(1, GLMBLen);
            wSum = 0;
            for idx = 1:GLMBLen
                newGLMB{idx} = GLMB{idx};
                wSum = wSum + newGLMB{idx}.w;
            end
            for idx = 1:GLMBLen
                newGLMB{idx}.w = newGLMB{idx}.w / wSum;
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function checkData(this, GLMB)
            wSum = 0;
            GLMBLen = length(GLMB);
            for idx = 1:GLMBLen
                wSum = wSum + GLMB{idx}.w;
            end
            if wSum ~= 1
                disp('Data error!');
                wSum
            end
        end
    end
end

