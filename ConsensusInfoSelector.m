classdef ConsensusInfoSelector < DistributedFusionAlgorithm
    %
    
    properties
        SensorNode;
    end
    
    methods
        % Constructor
        function td = ConsensusInfoSelector(SensorNode, L, theta)
            td.SensorNode = SensorNode;
            td.L = L;
            td.Time = 1;
            td.THETA = theta;
        end
        % Initialization
        function init(this)
            this.Time = 1;
            AgentNum = length(this.SensorNode);
            for idx = 1:AgentNum
                this.SensorNode{idx}.init();
            end
        end
        % Single step filtering and fusion
        function [TrackDistributions] = filter_update(this, Measure)
            AgentNum = length(this.SensorNode);
            % Prediction
            preBernoulliSets = cell(1, AgentNum);
            parfor idx = 1:AgentNum
                preBernoulliSets{idx} = this.SensorNode{idx}.timeUpdate(this.Time);
            end
            % Filtering
            estGLMB = cell(1, AgentNum);
            TotalNewbornSet = cell(1, AgentNum);
            TotalNewbornWeights = cell(1, AgentNum);
            parfor idx = 1:AgentNum
                [estGLMB{idx}, TotalNewbornSet{idx}, TotalNewbornWeights{idx}] =...
                    this.SensorNode{idx}.measurementUpdate(preBernoulliSets{idx}, Measure{idx});
            end
            for idx = 1:AgentNum
                this.SensorNode{idx}.CurrentGLMB = estGLMB{idx};
                this.SensorNode{idx}.NewbornWeights = TotalNewbornWeights{idx};
                this.SensorNode{idx}.CurrentNewborn = TotalNewbornSet{idx};
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             testLMB = this.SensorNode{1}.GLMB2LMB(estGLMB{1});
%             figure(2);
%             cla;
%             hold on
%             grid on
%             for idx = 1:length(preBernoulliSets{1})
%                 plot(preBernoulliSets{1}{idx}.x(1), preBernoulliSets{1}{idx}.x(3), 'Marker', '+');
%                 tmpLabelTxt = sprintf('(%d,%d)', preBernoulliSets{1}{idx}.l(1), preBernoulliSets{1}{idx}.l(2));
%                 text(preBernoulliSets{1}{idx}.x(1), preBernoulliSets{1}{idx}.x(3), tmpLabelTxt);
%             end
%             for idx = 1:length(testLMB)
%                 plot(testLMB{idx}.x(1), testLMB{idx}.x(3),'Marker','*');
%                 tmpLabelTxt = sprintf('(%d,%d)', testLMB{idx}.l(1),testLMB{idx}.l(2));
%                 text(testLMB{idx}.x(1), testLMB{idx}.x(3), tmpLabelTxt);
%             end
%             hold off
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             pre_r_sum = zeros(1, AgentNum);
%             for idx1 = 1:AgentNum
%                 tmp_r_sum = 0;
%                 preBernoulliSet = preBernoulliSets{idx1};
%                 for idx2 = 1:length(preBernoulliSet)
%                     tmp_r_sum = tmp_r_sum + preBernoulliSet{idx2}.r;
%                 end
%                 pre_r_sum(idx1) = tmp_r_sum;
%             end
%             pre_r_sum
%             ave_r_sum = 0;
%             for idx1 = 1:AgentNum
%                 tmp_r_sum = 0;
%                 preBernoulliSet = preBernoulliSets{idx1};
%                 for idx2 = 1:length(preBernoulliSet)
%                     tmp_r_sum = tmp_r_sum + preBernoulliSet{idx2}.r;
%                 end
%                 ave_r_sum = ave_r_sum + tmp_r_sum;
%             end
%             ave_r_sum = ave_r_sum / AgentNum;
%             ave_r_sum
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %%%%%%%%%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             disp('Filtering');
%             N = zeros(1, AgentNum);
%             for agt = 1:AgentNum
%                 tmpGLMBset = estGLMB{agt};
%                 tempILen = length(tmpGLMBset);
%                 tmpN = 0;
%                 for i = 1:tempILen
%                     tmpN = tmpN +  tmpGLMBset{i}.w * size(tmpGLMBset{i}.Label, 1);
%                 end
%                 N(agt) = tmpN;
%             end
%             N
%             sum(N)
%             disp('Local Newborn Density:');
%             NewbornR = zeros(1, AgentNum);
%             for idx = 1:AgentNum
%                 tmpR = 0;
%                 tmpLen = length(TotalNewbornSet{idx});
%                 for i = 1:tmpLen
%                     tmpR = tmpR + TotalNewbornSet{idx}{i}.r;
%                 end
%                 NewbornR(idx) = tmpR;
%             end
%             NewbornR
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fusion
            for l = 1:this.L
                % Obtain the output data of each agent
                TotalTransportDataGLMB = cell(1, AgentNum);
                TotalTransportDataNewborn = cell(1, AgentNum);
                TotalNeighbour = cell(1, AgentNum);
                parfor i = 1:AgentNum
                    TotalTransportDataGLMB{i} = this.SensorNode{i}.genTransportDataGLMB();
                    TotalTransportDataNewborn{i} = this.SensorNode{i}.genTransportDataNewborn();
                    TotalNeighbour{i} = this.SensorNode{i}.getNeighbour();
                end
                % Obtain the received data for each agent
                TotalGLMBsets = cell(1, AgentNum);
                TotalNewbornSets = cell(1, AgentNum);
                for i = 1:AgentNum
                    % i: the index of the received agent
                    NeighbourLen = length(TotalNeighbour{i});
                    estGLMB = cell(1, NeighbourLen);
                    NewbornSets = cell(1, NeighbourLen);
                    for idx1 = 1:NeighbourLen
                        tmpIndex = TotalNeighbour{i}(idx1);
                        estGLMB{idx1} = TotalTransportDataGLMB{tmpIndex};
                        NewbornSets{idx1} = TotalTransportDataNewborn{tmpIndex};
                    end
                    TotalGLMBsets{i} = estGLMB;
                    TotalNewbornSets{i} = NewbornSets;
                end
                % Local fusion
                CurrentFusedGLMB = cell(1, AgentNum);
                CurrentUpdateNewborn = cell(1, AgentNum);
                parfor i = 1:AgentNum
                    CurrentFusedGLMB{i} = this.SensorNode{i}.fusion(TotalGLMBsets{i});
                    CurrentUpdateNewborn{i} = this.SensorNode{i}.concateNewbornSet(TotalNewbornSets{i});
                end
                for idx1 = 1:AgentNum
                    this.SensorNode{idx1}.CurrentGLMB = CurrentFusedGLMB{idx1};
                    this.SensorNode{idx1}.CurrentNewborn = CurrentUpdateNewborn{idx1};
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             TotalWeightsSum = zeros(1, AgentNum);
%             for agt = 1:AgentNum
%                 weights = zeros(1, length(this.SensorNode{agt}.CurrentGLMB));
%                 for idx = 1:length(this.SensorNode{agt}.CurrentGLMB)
%                     weights(idx) = this.SensorNode{agt}.CurrentGLMB{idx}.w;
%                 end
%                 TotalWeightsSum(agt) = sum(weights);
%             end
%             TotalWeightsSum
%             disp('Newborn density after combining:');
%             TotalRSum = zeros(1, AgentNum);
%             TotalRLen = zeros(1, AgentNum);
%             for agt = 1:AgentNum
%                 tmpNewbornSet = this.SensorNode{agt}.NewbornSet;
%                 TotalRLen(agt) = length(tmpNewbornSet);
%                 tmp_r_sum = 0;
%                 for idx = 1:length(tmpNewbornSet)
%                     tmp_r_sum = tmp_r_sum + tmpNewbornSet{idx}.r;
%                 end
%                 TotalRSum(agt) = tmp_r_sum;
%             end
% %             ave_r_sum = sum(TotalRSum) / length(TotalRSum);
%             TotalRSum
% %             ave_r_sum
%             ave_r_sum = 0;
%             for agt = 1:AgentNum
%                 tmpNewbornSet = this.SensorNode{agt}.NewbornSet;
%                 tmp_r_sum = 0;
%                 for idx = 1:length(tmpNewbornSet)
%                     tmp_r_sum = tmp_r_sum + tmpNewbornSet{idx}.r;
%                 end
%                 ave_r_sum = ave_r_sum + tmp_r_sum;
%             end
%             ave_r_sum = ave_r_sum / AgentNum;
%             TotalRLen
%             figure(2);
%             hold on
%             testNewbornSet = this.SensorNode{1}.NewbornSet;
%             testNewbornSetLen = length(testNewbornSet);
%             for idx = 1:testNewbornSetLen
%                 plot(testNewbornSet{idx}.x(1), testNewbornSet{idx}.x(3),'Marker','o','Color','b');
%             end
%             hold off
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Synthesizing data
            TrackDistributions = cell(1, AgentNum);
            UpdateNewbornSets = cell(1, AgentNum);
            parfor idx = 1:AgentNum
                TrackDistributions{idx} = this.SensorNode{idx}.GLMB2LMB(this.SensorNode{idx}.CurrentGLMB);
                UpdateNewbornSets{idx} = this.SensorNode{idx}.mergingNewbornSet2();
            end
            for idx = 1:AgentNum
                this.SensorNode{idx}.NewbornWeights = UpdateNewbornSets{idx};
                this.SensorNode{idx}.BernoulliSet = TrackDistributions{idx};
            end
            %%%%%%%%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             disp('Fusion');
%             rSum = 0;
%             testNewbornSet = UpdateNewbornSets{1};
%             testNewbornSetLen = length(testNewbornSet);
%             testNewbornSetLen
%             for idx = 1:testNewbornSetLen
%                 rSum = rSum + testNewbornSet{idx}.r;
%             end
%             rSum
%             rSumFused = 0;
%             testBernoulliSet = TrackDistributions{1};
%             testBernoulliSetLen = length(testBernoulliSet);
%             for idx = 1:testBernoulliSetLen
%                 rSumFused = rSumFused + testBernoulliSet{idx}.r;
%             end
%             rSumFused
%             N = 0;
%             testBernoulliSet = TrackDistributions{1};
%             testBernoulliSetLen = length(testBernoulliSet);
%             for idx = 1:testBernoulliSetLen
%                 if testBernoulliSet{idx}.r > 0.5
%                     N = N + 1;
%                 end
%             end
%             N
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        % Filtering
        function [TotalTrack, ER, EN] = filter(this, measure)
            AgentNum = length(this.SensorNode);
            TimeLength = length(measure{1});
            TotalTrack = cell(1, AgentNum);
            if AgentNum ~= length(measure)
                disp('The number of agents is not valid!!!');
            end
            this.init();
            TotalTrackDistributions = cell(1, TimeLength);
            ER = zeros(TimeLength, 1);
            EN = zeros(TimeLength, 1);
            for k = 1:TimeLength
%                 str = sprintf('k=%d',k);
%                 disp(str);
                tempMeasure = cell(1, AgentNum);
                parfor i = 1:AgentNum
                    tempMeasure{i} = measure{i}{k};
                end
                TrackDistributions = this.filter_update(tempMeasure);
                ER(k) = this.getN2(TrackDistributions{1});
                EN(k) = this.getN(TrackDistributions{1});
                TotalTrackDistributions{k} = TrackDistributions;
                this.Time = this.Time + 1;
            end
            % Obtain the tracks for each agent
            parfor i = 1:AgentNum
                Track = cell(0,0);
                time = 1;
                for k = 1:TimeLength
                    Track = this.updateTrack(Track, TotalTrackDistributions{k}{i}, time);
                    time = time + 1;
                end
                TotalTrack{i} = Track;
            end
        end
        function TimeUpdate(this)
            this.Time = this.Time + 1;
        end
    end
    
    methods
        % Update the track set
        function fusedTrack = updateTrack(this, oldTrack, fusedDistribution, time)
            fusedDistributionLen = length(fusedDistribution);
            r = zeros(1, fusedDistributionLen);
            parfor k = 1:fusedDistributionLen
                r(k) = fusedDistribution{k}.r;
            end
            estimateNumber = length(find(r >= 0.5));
            estimateNumber = min(estimateNumber, fusedDistributionLen);
            [~, index] = sort(r, 'descend');
            reorderDistribution = cell(1, fusedDistributionLen);
            parfor k = 1:fusedDistributionLen
                reorderDistribution{k} = fusedDistribution{index(k)};
            end
            TrackLen = length(oldTrack);
            fusedTrack = oldTrack;
            count = TrackLen + 1;
            for i = 1:estimateNumber
                flag = true;
                for j = 1:TrackLen
                    if isequal(oldTrack{j}.l, reorderDistribution{i}.l)
                        tempX = reorderDistribution{i}.x;
                        fusedTrack{j}.X = [oldTrack{j}.X, tempX];
                        fusedTrack{j}.t = [oldTrack{j}.t, time];
                        flag = false;
                    end
                end
                if flag
                    tempX = reorderDistribution{i}.x;
                    fusedTrack{count}.X = tempX;
                    fusedTrack{count}.t = time;
                    fusedTrack{count}.l = reorderDistribution{i}.l;
                    count = count + 1;
                end
            end
        end
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
        function N = getN2(this, LMB)
            LMBLen = length(LMB);
            N = zeros(1, LMBLen);
            parfor j = 1:LMBLen
                N(j) = LMB{j}.r;
            end
            N = sum(N);
        end
    end
end

