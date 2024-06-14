classdef CMIL_GLMB_DF < FusionAlgorithm
    % The implementation of the fusion algorithm proposed in Reference:
    % Reference:
    % Lin Gao, Giorgio Battistelli, Luigi Chisci, "Fusion of Labeled RFS
    % Densities with Different Fields of View", IEEE Transactions on
    % aerospace and electronic systems, vol. 58, no. 6, pp. 5908-5924,
    % 2022.
    % Fusion structure part
    
    properties
        Time;
        AdaptiveWeightFlag;
        L;
    end
    
    methods
        % Constructor
        function [td] = CMIL_GLMB_DF(fusionWeights, filter, adaptiveWeightFlag, l)
            td.FusionWeights = fusionWeights;
            td.SensorNode = filter;
            td.Time = 1;
            td.AdaptiveWeightFlag = adaptiveWeightFlag;
            td.L = l;
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
            for idx = 1:AgentNum
                this.SensorNode{idx}.IndexCount = length(this.SensorNode{idx}.NewbornSet) + 1;
            end
            % Filtering
            estGLMB = cell(1, AgentNum);
            TotalNewbornSet = cell(1, AgentNum);
            parfor idx = 1:AgentNum
                [estGLMB{idx}, TotalNewbornSet{idx}] = this.SensorNode{idx}.measurementUpdate(preBernoulliSets{idx}, Measure{idx});
            end
            for idx = 1:AgentNum
                this.SensorNode{idx}.CurrentGLMB = estGLMB{idx};
                this.SensorNode{idx}.NewbornSet = TotalNewbornSet{idx};
            end
            % Fusion
            for l = 1:this.L
                % Obtain the output data of each agent
                TotalTransportDataGLMB = cell(1, AgentNum);
                TotalNeighbour = cell(1, AgentNum);
                parfor i = 1:AgentNum
                    TotalTransportDataGLMB{i} = this.SensorNode{i}.genTransportDataGLMB();
                    TotalNeighbour{i} = this.SensorNode{i}.getNeighbour();
                end
                % Obtain the received data for each agent
                TotalGLMBsets = cell(1, AgentNum);
                for i = 1:AgentNum
                    % i: the index of the received agent
                    NeighbourLen = length(TotalNeighbour{i});
                    estGLMB = cell(1, NeighbourLen);
                    for idx1 = 1:NeighbourLen
                        tmpIndex = TotalNeighbour{i}(idx1);
                        estGLMB{idx1} = TotalTransportDataGLMB{tmpIndex};
                    end
                    TotalGLMBsets{i} = estGLMB;
                end
                % Local fusion
                CurrentFusedGLMB = cell(1, AgentNum);
                parfor i = 1:AgentNum
                    CurrentFusedGLMB{i} = this.SensorNode{i}.fusion(TotalGLMBsets{i}, this.AdaptiveWeightFlag, this.Time);
                end
                for idx1 = 1:AgentNum
                    this.SensorNode{idx1}.CurrentGLMB = CurrentFusedGLMB{idx1};
                end
            end
            % Synthesizing data
            TrackDistributions = cell(1, AgentNum);
            parfor idx = 1:AgentNum
                tmpDistribution = this.SensorNode{idx}.GLMB2LMB(this.SensorNode{idx}.CurrentGLMB);
                TrackDistributions{idx} = this.SensorNode{idx}.trackPruning(tmpDistribution);
            end
            for idx = 1:AgentNum
                this.SensorNode{idx}.BernoulliSet = TrackDistributions{idx};
            end
            %%%%%%%%%%%%%%%%%%%%%%%%% TEST %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             disp('Fusion');
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

