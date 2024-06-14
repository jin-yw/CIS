classdef Sensor < handle
    % Sensor model
    
    properties
        MeasurementModel;
        KinematicModel;
        Pd;
        Lambda;
        ClutterModel;
        Range;
        SensorID;
        Neighbour;
        BernoulliSet;
        NewbornWeights;
        NewbornDistribution;
        Ps;
        Pg;
        Gamma;
        THETA;
        N;
        CurrentGLMB;
        CurrentNewborn;
        rBMax;
        LambdaB;
        U_star;
    end
    
    methods
        % Constructor
        function td = Sensor(measurementModel, kinematicModel, clutterModel, pd, ps, pg, gamma, theta, lambda, range, sensorID, neighbour, newbornDistribution, n, rbmax, lambdab)
            td.MeasurementModel = measurementModel;
            td.KinematicModel = kinematicModel;
            td.Pd = pd;
            td.Lambda = lambda;
            td.ClutterModel = clutterModel;
            td.Range = range;
            td.SensorID = sensorID;
            td.Neighbour = neighbour;
            td.BernoulliSet = cell(0,0);
            td.NewbornDistribution = newbornDistribution;
            [Nx, Ny] = size(newbornDistribution);
            td.NewbornWeights = zeros(1, Nx*Ny);
            td.Ps = ps;
            td.Pg = pg;
            td.Gamma = gamma;
            td.THETA = theta;
            td.N = n;
            td.CurrentGLMB = cell(0,0);
            td.CurrentNewborn = cell(0,0);
            td.rBMax = rbmax;
            td.LambdaB = lambdab;
            u_star = zeros(Nx*Ny, Nx*Ny);
            for idx1 = 1:(Nx*Ny)
                for idx2 = 1:(Nx*Ny)
                    j1 = mod(idx1, Ny);
                    i1 = floor(idx1/Ny) + 1;
                    if j1 == 0
                        j1 = Ny;
                        i1 = i1-1;
                    end
                    i2 = floor(idx2/Ny) + 1;
                    j2 = mod(idx2, Ny);
                    if j2 == 0
                        j2 = Ny;
                        i2 = i2 - 1;
                    end
                    m1 = newbornDistribution{i1, j1}.x;
                    P1 = newbornDistribution{i1, j1}.P;
                    try
                        m2 = newbornDistribution{i2, j2}.x;
                    catch
                        disp('Stop here!');
                    end
                    P2 = newbornDistribution{i2, j2}.P;
                    P = P1 + P2;
                    uu = 1/sqrt(det(2*pi*P)) * exp(-1/2 * (m1-m2)' * pinv(P) * (m1-m2));
                    u_star(idx1, idx2) = uu;
                end
            end
            td.U_star = u_star;
            td.setPath();
        end
        function setPath(this)
            addpath('_common\');
        end
        % Initialization
        function init(this)
            this.NewbornWeights = zeros(size(this.NewbornDistribution));
            this.BernoulliSet = cell(0,0);
            this.CurrentGLMB = cell(0,0);
            this.CurrentNewborn = cell(0,0);
        end
        % Time update
        function [preBernoulliSet] = timeUpdate(this, Time)
            % Prediction
            oldBernoulliLen = length(this.BernoulliSet);
            preBernoulliSetOld = cell(1, oldBernoulliLen);
            % prediction for existed targets
            parfor idx = 1:oldBernoulliLen
                preBernoulliSetOld{idx} = this.updateBernoulliComponent(this.BernoulliSet{idx});
            end
            % prediction for newborn targets
            preBernoulliSetNew = this.genNewbornTarget(Time);
            preBernoulliSet = this.combineBernoulliSet(preBernoulliSetOld, preBernoulliSetNew);
        end
        % Measurement update
        function [GLMB, newbornSet, newbornWeights] = measurementUpdate(this, preBernoulliSet, measure)
            % Parallel Group Update
            preGLMB = this.LMB2GLMB(preBernoulliSet, this.N);
            [GLMB, newbornSet, newbornWeights] = this.getUpdateGLMB(preGLMB, measure);
            this.NewbornWeights = newbornWeights;
            this.CurrentNewborn = newbornSet;
            this.CurrentGLMB = GLMB;
        end
        function [fusedGLMB] = fusion(this, GLMBset)
            AgentNum = length(GLMBset);         % GLMBset{agent}{GLMB component}
            GLMBgroup = cell(0,0);
            count = 1;
            % Grouping
            for idx1 = 1:AgentNum
                tempGLMB = GLMBset{idx1};
                tempGLMBLen = length(tempGLMB);
                for idx2 = 1:tempGLMBLen
                    GLMBComponent = tempGLMB{idx2};
                    % GLMBComponent: Distribution, Label, w
                    tmpLabelSet = GLMBComponent.Label;
                    [index, reorderI] = this.findLabelSet(tmpLabelSet, GLMBgroup);
                    if index > 0
                        % Find the corresponding GLMB component set
                        wLen = length(GLMBgroup{index}.w);
                        GLMBgroup{index}.w(wLen + 1) = tempGLMB{idx2}.w;
                        if ~isempty(reorderI)
                            GLMBgroup{index}.w(wLen + 1) = tempGLMB{idx2}.w;
                            for idx3 = 1:length(reorderI)
                                GLMBgroup{index}.DistributionSet{wLen + 1}{idx3} = tempGLMB{idx2}.Distribution{reorderI(idx3)};
                            end
                        else
                            GLMBgroup{index}.DistributionSet{wLen + 1} = cell(0,0);
                        end
                    else
                        GLMBgroup{count}.LabelSet = tempGLMB{idx2}.Label;
                        GLMBgroup{count}.w = [tempGLMB{idx2}.w];
                        GLMBgroup{count}.DistributionSet = cell(0,0);
                        GLMBgroup{count}.DistributionSet{1} = tempGLMB{idx2}.Distribution;
                        count = count + 1;
                    end
                end
            end
            % Fusion
            GLMBgroupLen = length(GLMBgroup);
            fusedGLMB = cell(1, GLMBgroupLen);
            parfor idx = 1:GLMBgroupLen
                fusedGLMB{idx} = this.genFusedDistribution(GLMBgroup{idx});
            end
        end
        function updateNewbornSet = concateNewbornSet(this, TotalNewbornSet)
            AgentNum = length(TotalNewbornSet);
            GLMBgroup = cell(0,0);
            count = 1;
            % Grouping
            for idx1 = 1:AgentNum
                tempGLMB = TotalNewbornSet{idx1};
                tempGLMBLen = length(tempGLMB);
                for idx2 = 1:tempGLMBLen
                    GLMBComponent = tempGLMB{idx2};
                    % GLMBComponent: w, LabelSet
                    tmpLabelSet = GLMBComponent.LabelSet;
                    [index] = this.findLabelSet2(tmpLabelSet,GLMBgroup);
                    if index > 0
                        % Find the corresponding GLMB component set
                        wLen = length(GLMBgroup{index}.w);
                        GLMBgroup{index}.w(wLen + 1) = tempGLMB{idx2}.w;
                    else
                        GLMBgroup{count}.LabelSet = tempGLMB{idx2}.LabelSet;
                        GLMBgroup{count}.w = [tempGLMB{idx2}.w];
                        count = count + 1;
                    end
                end
            end
            % Fusion
            GLMBgroupLen = length(GLMBgroup);
            updateNewbornSet = cell(1, GLMBgroupLen);
            parfor idx = 1:GLMBgroupLen
                maxW = max(GLMBgroup{idx}.w);
                updateNewbornSet{idx}.w = maxW;
                updateNewbornSet{idx}.LabelSet = GLMBgroup{idx}.LabelSet;
            end
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
                if r > 0
                    x = x / r;
                    for i = 1:bernoulliSetLen
                        P = P + bernoulliSet{i}.r * (bernoulliSet{i}.P + (bernoulliSet{i}.x - x) * (bernoulliSet{i}.x - x)');
                    end
                    P = P / r;
                    LMB{idx}.x = x;
                    LMB{idx}.P = P;
                    LMB{idx}.l = LMBgroup{idx}.l;
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
            LMB = this.trackPruning(LMB);
            this.BernoulliSet = LMB;
        end
        function [UpdateNewbornWeights] = mergingNewbornSet2(this)
            % GLMB->LMB
            GLMBLen = length(this.CurrentNewborn);
            [Nx,Ny] = size(this.NewbornDistribution);
            Len = Nx * Ny;
            UpdateNewbornWeights = zeros(1, Len);
            % Normalization
            wSum = 0;
            for idx = 1:GLMBLen
                wSum = wSum + this.CurrentNewborn{idx}.w;
            end
            for idx = 1:GLMBLen
                this.CurrentNewborn{idx}.w = this.CurrentNewborn{idx}.w / wSum;
            end
            % Marginalizing
            for idx = 1:GLMBLen
                LabelSet = this.CurrentNewborn{idx}.LabelSet;
                LabelSetLen = length(LabelSet);
                for i = 1:LabelSetLen
                    UpdateNewbornWeights(LabelSet(i)) = UpdateNewbornWeights(LabelSet(i)) + this.CurrentNewborn{idx}.w;
                end
            end
            % 
            this.NewbornWeights = UpdateNewbornWeights;
        end
        function [transportDataGLMB] = genTransportDataGLMB(this)
            transportDataGLMB = this.CurrentGLMB;
        end
        function [newbornSet] = genTransportDataNewborn(this)
            newbornSet = this.CurrentNewborn;
        end
        function neighbour = getNeighbour(this)
            neighbour = this.Neighbour;
        end
    end
    
    methods
        % Prediction for existed targets
        function bernoulliComponent = updateBernoulliComponent(this, oldBernoulliComponent)
            bernoulliComponent.r = this.Ps * oldBernoulliComponent.r;
            bernoulliComponent.l = oldBernoulliComponent.l;
            bernoulliComponent.x = this.KinematicModel.estimateNewState(oldBernoulliComponent.x);
            F = this.KinematicModel.getLinear(oldBernoulliComponent.x);
            bernoulliComponent.P = F * oldBernoulliComponent.P * F' + this.KinematicModel.CovMatrix;
        end
        % prediction for newborn targets
        function [newbornBernoulliSets] = genNewbornTarget(this, time)
            Len = length(this.NewbornWeights);
            [Nx, Ny] = size(this.NewbornDistribution);
            newbornBernoulliSets = cell(0,0);
            count = 1;
            for idx1 = 1:Nx
                for idx2 = 1:Ny
                    index = (idx1-1) * Ny + idx2;
                    bernoulliComponent.x = this.NewbornDistribution{idx1, idx2}.x;
                    bernoulliComponent.P = this.NewbornDistribution{idx1, idx2}.P;
                    bernoulliComponent.l = [time, index];
                    bernoulliComponent.r = this.NewbornWeights(index);
                    if bernoulliComponent.r > 0.001
                        newbornBernoulliSets{count} = bernoulliComponent;
                        count = count + 1;
                    end
                end
            end
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
        function [GLMB, newbornSet, newbornWeights] = getUpdateGLMB(this, preGLMB, measure)
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
                ETA = [ETA; tempETA];       % 1: Measurement Number; 2: Distribution Number
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
            [newbornSet, newbornWeights] = this.genBirthDistribution(posteriorGLMB, measure);
            GLMB = this.simplifiedGLMB(posteriorGLMB, Distribution, Label, Z);
        end
        % Single target update 2
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
        function [newbornSet, newbornWeights] = genBirthDistribution(this, GLMB, measure)
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
                newbornSet{idx}.P = diag([50; 50; 50; 50].^2);
                newbornSet{idx}.r = min(this.rBMax, rB(idx) * this.LambdaB);
                newbornSet{idx}.r = double(vpa(newbornSet{idx}.r, 7));
            end
            % 剪枝
            tempNewbornSet = cell(0,0);
            count = 1;
            for idx = 1:measureLen
                if newbornSet{idx}.r > this.THETA
                    tempNewbornSet{count} = newbornSet{idx};
                    count = count + 1;
                end
            end
            % Approximation
            [Nx, Ny] = size(this.NewbornDistribution);
            Len1 = Nx * Ny;
            Len2 = length(tempNewbornSet);
            V = zeros(Len2, Len1);
            for idx1 = 1:Len2
                for idx2 = 1:Len1
                    x1 = tempNewbornSet{idx1}.x;
                    P1 = tempNewbornSet{idx1}.P;
                    i2 = floor(idx2/Ny) + 1;
                    j2 = mod(idx2, Ny);
                    if j2 == 0
                        j2 = Ny;
                        i2 = i2 - 1;
                    end
                    x2 = this.NewbornDistribution{i2, j2}.x;
                    P2 = this.NewbornDistribution{i2, j2}.P;
                    P = P1 + P2;
                    vv = 1/sqrt(det(2*pi*P)) * exp(-1/2 * (x1-x2)' * pinv(P) * (x1-x2));
                    V(idx1, idx2) = vv;
                end
            end
            r = zeros(1, Len2);
            for idx = 1:Len2
                r(idx) = tempNewbornSet{idx}.r;
            end
            f = -r * V;
            opts = optimset('Algorithm','active-set','Display','off');
            x0 = ones(Len1,1) * this.rBMax;
            tmpMatrix = [this.U_star; f];
            tmpAbsMatrix = abs(tmpMatrix);
            tmpVal = mean(tmpAbsMatrix(:));
            tmpVal = 1 / tmpVal * 100;
            [x, fval] = quadprog(this.U_star * tmpVal, f * tmpVal, [], [], [], [], ones(Len1,1) * 0, ones(Len1,1) * 1, x0, opts);
            newbornWeights = x';
            newX = [];
            newIDX = [];
            for idx = 1:Len1
                if x(idx) > this.THETA
                    newX = [newX, x(idx)];
                    newIDX = [newIDX, idx];
                end
            end
            x = newX;
            Len3 = length(x);
            % LMB->GLMB
            C1 = zeros(Len3, Len3);
            C2 = zeros(Len3, Len3);
            for idx1 = 1:Len3
                for idx2 = 1:Len3
                    if idx1 == idx2
                        C1(idx1, idx2) = -log(x(idx1));
                        C2(idx1, idx2) = -log(1 - x(idx1));
                    else
                        C1(idx1, idx2) = inf;
                        C2(idx1, idx2) = inf;
                    end
                end
            end
            C = [C1, C2];
            [assignments, cost] = murty_custom(C, 30);
            tmpLen = length(cost);
            newbornSet = cell(1, tmpLen);
            sumW = 0;
            for idx = 1:tmpLen
                tmpW = exp(-cost(idx));
                newbornSet{idx}.w = tmpW;
                sumW = sumW + tmpW;
                tmpLabelSet = [];
                for i = 1:Len3
                    tmpIndex = assignments(idx, i);
                    if tmpIndex < Len3
                        tmpLabelSet = [tmpLabelSet, newIDX(i)];
                    end
%                     tmpVec = xSet{idx}(i,:);
%                     tmpIndex = find(tmpVec == 1);
%                     if tmpIndex(1) < Len3
%                         tmpLabelSet = [tmpLabelSet, newIDX(i)];
%                     end
                end
                newbornSet{idx}.LabelSet = tmpLabelSet;
            end
            % Normalization
            for idx = 1:tmpLen
                newbornSet{idx}.w = newbornSet{idx}.w / sumW;
            end
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
                        GLMB{idx}.Z{k} = Z(tempI(k));
                    end
                end
            end
            GLMBLen = length(GLMB);
            for idx = 1:GLMBLen
                GLMB{idx}.w = GLMB{idx}.w / wSum;
                GLMB{idx}.w = double(vpa(GLMB{idx}.w, 7));
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
                LMB{idx} = oldLMB{index(idx)};
                if LMB{idx}.r >= 1
                    LMB{idx}.r = 1-eps;
                end
%                 LMB{idx}.r = vpa(LMB{idx}.r, 7);
            end
        end
        %%%%%%%%%%%%%%%%%%% FUNCTIONS FOR FUSION %%%%%%%%%%%%%%%%%%%%%%%%%%
        function [index, reorderI] = findLabelSet(this, refLabelSet, GLMBgroup)
            GLMBgroupLen = length(GLMBgroup);
            refLabelSetLen = size(refLabelSet,1);
            index = 0;
            reorderI = [];
            for idx = 1:GLMBgroupLen
                reorderI = [];
                LabelSet = GLMBgroup{idx}.LabelSet;
                LabelSetLen = size(LabelSet,1);
                if LabelSetLen ~= refLabelSetLen
                    continue;
                elseif (LabelSetLen == refLabelSetLen) && (refLabelSetLen == 0)
                    index = idx;
                else
                    flag = true;
                    tempIndex = 0;
                    for i = 1:LabelSetLen
                        tempFlag = false;
                        for j = 1:refLabelSetLen
                            if isequal(LabelSet(i,:), refLabelSet(j,:))
                                tempFlag = true;
                                tempIndex = j;
                                break;
                            end
                        end
                        if ~tempFlag
                            flag = false;
                            break;
                        end
                        reorderI = [reorderI, tempIndex];
                    end
                    if flag
                        index = idx;
                        break;
                    end
                end
            end
        end
        % The labels are represented by scalars.
        function [index] = findLabelSet2(this, refLabelSet, GLMBgroup)
            GLMBgroupLen = length(GLMBgroup);
            index = 0;
            for idx = 1:GLMBgroupLen
                LabelSet = GLMBgroup{idx}.LabelSet;
                tmpSet1 = intersect(refLabelSet, LabelSet);
                tmpSet2 = union(refLabelSet, LabelSet);
                if length(tmpSet1) == length(tmpSet2)
                    index = idx;
                    break;
                end
            end
        end
        function [fusedGLMBcomponent] = genFusedDistribution(this, GLMBComponentGroup)
            weights = GLMBComponentGroup.w;
            index = find(weights == max(weights));
            fusedGLMBcomponent.Distribution = GLMBComponentGroup.DistributionSet{index(1)};
            fusedGLMBcomponent.Label = GLMBComponentGroup.LabelSet;
            fusedGLMBcomponent.w = GLMBComponentGroup.w(index(1));
        end
    end
end

