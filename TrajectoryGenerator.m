classdef TrajectoryGenerator < handle
    % The base class of real data generator.
    
    properties
        InitPos;                % Initial positions of trajectories.
        ManeuverType;           % The type of maneuvers of objects.
        MeasurementType;        % The type of measurements.
        StartTime;              % The time that object appears.
        EndTime;                % The time that object ends.
        SimulationTime;         % The length of time of the simulation.
        Step;                   % The step of the simulation
        Pd;                     % The detection rate.
        Lambda;                 % The average number of clutters.
        Range;                  % The observation range.
    end
    
    methods
        %% Constructor
        function td = TrajectoryGenerator(initpos, maneuvertype, measurementtype, starttime, endtime, simulationtime, step, pd, lambda, range)
            td.InitPos = initpos;
            td.ManeuverType = maneuvertype;
            td.MeasurementType = measurementtype;
            td.StartTime = starttime;
            td.EndTime = endtime;
            td.SimulationTime = simulationtime;
            td.Step = step;
            td.Pd = pd;
            td.Lambda = lambda;
            td.Range = range;
        end
        
        %% Generate the tracks
        function Truth = generateTruth(this)
            time = 0:this.Step:this.SimulationTime;
            Truth.K = length(time);                         % length of data/number of scans.
            Truth.X = cell(Truth.K, 1);                     % ground truth for states of targets
            Truth.N = zeros(Truth.K, 1);                    % ground truth for number of targets
            Truth.L = cell(Truth.K, 1);                     % ground truth for labels of targets (k,i)
            Truth.track_list = cell(Truth.K, 1);            % absolute index target identities
            Truth.total_tracks = length(this.ManeuverType); % total number of appearing tracks
            for targetNum = 1:Truth.total_tracks
                targetState = this.InitPos{targetNum};
                tempIndex = find(time == this.StartTime(targetNum));
                tempTime = 0:this.Step:this.EndTime(targetNum);
                for k = tempIndex : min(length(tempTime), length(time))
                    Truth.X{k} = [Truth.X{k}, targetState];
                    Truth.track_list{k} = [Truth.track_list{k}, targetNum];
                    Truth.N(k) = Truth.N(k) + 1;
                    targetState = this.ManeuverType{targetNum}.estimateNewState(targetState);
                end
            end
        end
        
        %% Generate measurements
        function Measure = generateMeasure(this, Truth)
            Measure.K = Truth.K;
            Measure.Z = cell(Truth.K,1);
            for k = 1:Truth.K
                if Truth.N(k) > 0
                    for n = 1:size(Truth.X{k},2)
                        x = Truth.X{k}(:,n);
                        pd = this.MeasurementType.getPd(x);
                        temp = rand();
                        if temp <= pd
                            Measure.Z{k} = [Measure.Z{k}, this.MeasurementType.MeasurementModel.predictNewState(x)];
                        end
                    end
                    Nc = poissrnd(this.Lambda);
                    C = [];
                    for n = 1:size(this.Range, 1)
                        C = [C; unifrnd(this.Range(n,1), this.Range(n,2), 1, Nc)];
                    end
                    Measure.Z{k} = [Measure.Z{k}, C];
                end
            end
        end
%         function Measure = generateMeasure(this, Truth)
%             Measure.K = Truth.K;
%             Measure.Z = cell(Truth.K, 1);
%             for k = 1:Truth.K
%                 if Truth.N(k) > 0
%                     idx = find(rand(Truth.N(k), 1) <= this.Pd);
%                     for n = 1:length(idx)
%                         Measure.Z{k} = [Measure.Z{k}, this.MeasurementType.predictNewState(Truth.X{k}(:,idx(n)))];
%                     end
%                 end
%                 Nc = poissrnd(this.Lambda);
%                 C = [];
%                 for n = 1:size(this.Range, 1)
%                     C = [C; unifrnd(this.Range(n,1), this.Range(n,2), 1, Nc)];
%                 end
%                 Measure.Z{k} = [Measure.Z{k}, C];
%             end
%         end
    end
    
end

