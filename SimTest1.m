clear all;
close all;
clc;

%% Single step simulation
%% Trajectory settings
TimeLength = 100;
Step = 1;
TargetNum = 12;
Pd = 0.98;
SensorNum = 6;
Lambda = ones(1,SensorNum) * 10;
Time = 0 : Step : TimeLength;
Range = cell(1,SensorNum);
count = 1;
Range{count} = [0,1500;-pi, 0];
count = count + 1;
Range{count} = [0,1500;-pi, 0];
count = count + 1;
Range{count} = [0, 1500;-pi, 0];
count = count + 1;
Range{count} = [0, 1500; 0, pi];
count = count + 1;
Range{count} = [0, 1500; 0, pi];
count = count + 1;
Range{count} = [0, 1500; 0, pi];
initPos = cell(TargetNum, 1);
startTime = zeros(12,1);
endTime = zeros(12,1);
initPos{1} = [0;0;0;-10];           startTime(1) = 1;       endTime(1) = 70;
initPos{2} = [400; -10; -600; 5];   startTime(2) = 1;       endTime(2) = TimeLength + Step;
initPos{3} = [-800; 20; -200; -5];  startTime(3) = 1;       endTime(3) = 70;
initPos{4} = [400; -7; -600; -4];   startTime(4) = 20;      endTime(4) = TimeLength + Step;
initPos{5} = [400; -2.5; -600; 10]; startTime(5) = 20;      endTime(5) = TimeLength + Step;
initPos{6} = [0; 7.5; 0; -5];       startTime(6) = 20;      endTime(6) = TimeLength + Step;
initPos{7} = [-800; 12; -200; 7];   startTime(7) = 40;      endTime(7) = TimeLength + Step;
initPos{8} = [-200; 15; 800; -10];  startTime(8) = 40;      endTime(8) = TimeLength + Step;
initPos{9} = [-800; 3; -200; 15];   startTime(9) = 60;      endTime(9) = TimeLength + Step;
initPos{10} = [-200; -3; 800; -15]; startTime(10) = 60;     endTime(10) = TimeLength + Step;
initPos{11} = [0; -20; 0; -15];     startTime(11) = 80;     endTime(11) = TimeLength + Step;
initPos{12} = [-200; 15; 800; -5];  startTime(12) = 80;     endTime(12) = TimeLength + Step;
maneuverType = cell(TargetNum, 1);
StateStd = [10,10];
KinematicModel = CVmodel(Step, StateStd);
for idx = 1:TargetNum
    maneuverType{idx} = KinematicModel;
end
Std = cell(1,SensorNum);
count = 1;
Std{count} = [10; 0.00175*2 ];
count = count + 1;
Std{count} = [20; 0.00175*1];
count = count + 1;
Std{count} = [5; 0.00175*0.5];
count = count + 1;
Std{count} = [15; 0.00175*1.5];
count = count + 1;
Std{count} = [10; 0.00175*2];
count = count + 1;
Std{count} = [15; 0.00175*1];
SensorPos = cell(1,SensorNum);
count = 1;
SensorPos{count} = [-1000;1000];            % 1
count = count + 1;
SensorPos{count} = [0;1000];
count = count + 1;
SensorPos{count} = [1000;1000];             % 2
count = count + 1;
SensorPos{count} = [-1000; -1000];          % 3
count = count + 1;
SensorPos{count} = [0; -1000];
count = count + 1;
SensorPos{count} = [1000; -1000];           % 4
MeasurementModel = cell(1,SensorNum);
for idx = 1:SensorNum
    MeasurementModel{idx} = SphereSensor(Step, Std{idx}, SensorPos{idx}, class(KinematicModel));
end
agent = cell(1,SensorNum);
for idx = 1:SensorNum
    agent{idx} = Agent(Range{idx}, MeasurementModel{idx}, Pd);
end
TrajectoryGene = cell(1,SensorNum);
for idx = 1:SensorNum
    TrajectoryGene{idx} = TrajectoryGenerator(initPos, maneuverType, agent{idx}, startTime, endTime, TimeLength, Step, Pd, Lambda(idx), Range{idx});
end
Truth = cell(1,SensorNum);
for idx = 1:SensorNum
    Truth{idx} = TrajectoryGene{idx}.generateTruth();
end
Measure = cell(1,SensorNum);
for idx = 1:SensorNum
    Measure{idx} = TrajectoryGene{idx}.generateMeasure(Truth{idx});
end
measure = cell(1,SensorNum);
for idx = 1:SensorNum
    measure{idx} = Measure{idx}.Z;
end
% Extrack trakcs
X_track = cell(TargetNum, 1);
for k = 1:Truth{1}.K
    tempTrackList = Truth{1}.track_list{k};
    for n = 1:length(tempTrackList)
        X_track{tempTrackList(n)} = [X_track{tempTrackList(n)}, [Time(k); Truth{1}.X{k}(:,n)]];
    end
end
figure(1);
hold on
grid on
TrackLen = length(X_track);
for idx = 1:TrackLen
    plot(X_track{idx}(2,:),X_track{idx}(4,:),'color','k', 'Marker', '.');
end
SensorRangeLine = cell(1,SensorNum);
for idx = 1:SensorNum
    tempRangeLine = zeros(2,2000);
    distance = Range{idx}(2,2) - Range{idx}(2,1);
    distance = distance / 2000;
    for k = 1:2000
        tempX = Range{idx}(1,2) * cos(Range{idx}(2,1) + distance * k);
        tempY = Range{idx}(1,2) * sin(Range{idx}(2,1) + distance * k);
        tempRangeLine(1,k) = tempX;
        tempRangeLine(2,k) = tempY;
    end
    tempRangeLine = [zeros(2,1), tempRangeLine, zeros(2,1)];
    SensorRangeLine{idx} = bsxfun(@plus, tempRangeLine, SensorPos{idx});
end
for idx = 1:SensorNum
    plot(SensorRangeLine{idx}(1,:), SensorRangeLine{idx}(2,:), 'color', 'k');
end

for idx = 1:SensorNum
    patch('XData', SensorRangeLine{idx}(1,:), 'YData', SensorRangeLine{idx}(2,:));
end
alpha(0.05);
hold off
% obtain the overlapping indicator value
TotalRangeX = [-10000, 10000; 0,0; -10000,10000;0,0];
AreaSize = 20000*20000;
overlappingRate = overlappingRateCalc(Range, agent, TotalRangeX, AreaSize);
overlappingRate
%% Filter and Fusion algorithms
Nx = 11;
Ny = 11;
NewbornDistribution = cell(Nx, Ny);
Dy = 2000 / (Nx - 1);
Dx = 2000 / (Ny - 1);
for idx1 = 1:Nx
    for idx2 = 1:Ny
        x = [-1000; 0; 1000; 0] + [(idx2-1) * Dx; 0; -(idx1-1) * Dy; 0];
        P = diag([Dx/6;Dx/6;Dx/6;Dx/6].^2);
        NewbornDistribution{idx1, idx2}.x = x;
        NewbornDistribution{idx1, idx2}.P = P;
    end
end
% Display the NewbornDistribution
% figure(10);
% SamplingPointsNum = 10000;
% Points = zeros(2, SamplingPointsNum);
% for idx = 1:SamplingPointsNum
%     xx = randi(Nx);
%     yy = randi(Ny);
%     pp = mvnrnd(NewbornDistribution{xx,yy}.x, NewbornDistribution{xx,yy}.P);
%     Points(:,idx) = pp(:,[1,3]);
% end
% scatter(Points(1,:), Points(2,:));
% hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ClutterModel = cell(1, SensorNum);
for idx = 1:SensorNum
    ClutterModel{idx} = GeneratorUniform(SensorPos{idx}, Range{idx});
end
Pb = 0.28;
Ps = 0.85;
Pg = 0.95;
Gamma = chi2inv(Pg, 2);
THETA = 1e-5;
N = 180;
LambdaB = 1;
rBMax = 0.2;
L2 = 3;
L1 = 1;
% Sensors
SensorNode1 = cell(1, SensorNum);
Neighbours = cell(1, SensorNum);
Neighbours{1} = [1,2,4];
Neighbours{2} = [1,2,3];
Neighbours{3} = [2,3,6];
Neighbours{4} = [1,4,5];
Neighbours{5} = [4,5,6];
Neighbours{6} = [3,5,6];
for idx = 1:SensorNum
    SensorNode1{idx} = Sensor(MeasurementModel{idx}, KinematicModel, ClutterModel{idx}, Pd, Ps, Pg, Gamma, THETA, Lambda(idx), Range{idx}, idx, Neighbours{idx}, NewbornDistribution, N, rBMax, LambdaB);
end
% CMIL-parameters
NeighboursCMIL = cell(1, SensorNum);
% NeighboursCMIL{1} = [2,3,4,5,6];
% NeighboursCMIL{2} = [1,3,4,5,6];
% NeighboursCMIL{3} = [1,2,4,5,6];
% NeighboursCMIL{4} = [1,2,3,5,6];
% NeighboursCMIL{5} = [1,2,3,4,6];
% NeighboursCMIL{6} = [1,2,3,4,5];
NeighboursCMIL{1} = [2,4];
NeighboursCMIL{2} = [1,3];
NeighboursCMIL{3} = [2,6];
NeighboursCMIL{4} = [1,5];
NeighboursCMIL{5} = [4,6];
NeighboursCMIL{6} = [3,5];
fusionWeights = cell(1, SensorNum);
for idx = 1:SensorNum
    tmpLen = length(NeighboursCMIL{idx}) + 1;
    fusionWeights{idx} = ones(1, tmpLen) / tmpLen;
end
SensorNodeCMIL = cell(1, SensorNum);
for idx = 1:SensorNum
    SensorNodeCMIL{idx} = SensorCMIL(MeasurementModel{idx}, KinematicModel, Pd, Lambda(idx), ClutterModel{idx}, Range{idx}, ...
        idx, NeighboursCMIL{idx}, Ps, Pg, Gamma, THETA, N, rBMax, LambdaB, fusionWeights{idx}, 30);
end
% IS
disp('IS');
FusionInfoSelector = InfoSelector(KinematicModel, agent, Ps, Pg, Gamma, Lambda, THETA, NewbornDistribution, ClutterModel, N, rBMax, LambdaB);
TracksInfoSelector = FusionInfoSelector.filter(measure);
% CIS1
disp('CIS');
FusionConsensusInfoSelector1 = ConsensusInfoSelector(SensorNode1, L1, THETA);
TotalTracksConsensusInfoSelector1 = FusionConsensusInfoSelector1.filter(measure);
TracksConsensusInfoSelector1 = TotalTracksConsensusInfoSelector1{1};
% CIS2
disp('CIS');
FusionConsensusInfoSelector2 = ConsensusInfoSelector(SensorNode1, L2, THETA);
TotalTracksConsensusInfoSelector2 = FusionConsensusInfoSelector2.filter(measure);
TracksConsensusInfoSelector2 = TotalTracksConsensusInfoSelector2{1};
% CMIL
FusionCMIL = CMIL_GLMB_DF(fusionWeights, SensorNodeCMIL, true, L2);
TotalTracksCMIL = FusionCMIL.filter(measure);
TracksCMIL = TotalTracksCMIL{1};
%% Plotting data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
colorBoard = {[1,0,1],...
    [1,0,0.5],...
    [0, 0.5, 1],...
    [0, 0.5, 0],...
    [0.4, 0.5, 0.8],...
    [0.4, 1, 0.8],...
    [0.7, 1, 0.7],...
    [1, 0.5, 0.5],...
    [1, 1, 0],...
    [0.6, 0.6, 0.6],...
    [0.7, 0.7, 0.7],...
    [0.5, 0, 0.5]};
% IS
figure(2);
hold on;
subplot(211); box on;

for i = 1:Truth{1}.total_tracks
    hline1 = line(X_track{i}(1,:),X_track{i}(2,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0 * ones(1,3));
end

for k = 1:length(TracksInfoSelector)
    colorIndex = mod(k,length(colorBoard));
    if colorIndex == 0
        colorIndex = length(colorBoard);
    end
    if size(TracksInfoSelector{k}.X,2) > 0
        tempTime = Time(TracksInfoSelector{k}.t);
        hline2 = line(tempTime, TracksInfoSelector{k}.X(1,:),'LineStyle','-','Marker','.','Markersize',8,'Color',colorBoard{colorIndex});
    end
end

subplot(212);box on;


for i = 1:Truth{1}.total_tracks
    yhline1 = line(X_track{i}(1,:), X_track{i}(4,:), 'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3));
end

for k = 1:length(TracksInfoSelector)
    colorIndex = mod(k,length(colorBoard));
    if colorIndex == 0
        colorIndex = length(colorBoard);
    end
    if size(TracksInfoSelector{k}.X,2) > 0
        tempTime = Time(TracksInfoSelector{k}.t);
        hline2 = line(tempTime, TracksInfoSelector{k}.X(3,:),'LineStyle','-','Marker','.','Markersize',8,'Color',colorBoard{colorIndex});
    end
end

subplot(211); xlabel('Time:s'); ylabel('x:m');
set(gca, 'XLim', [1 Truth{1}.K]); set(gca,'YLim',[-1000,1000]);
legend([hline2 hline1],'Estimation','Ground truth');

subplot(212); xlabel('Time:s'); ylabel('y:m');
set(gca, 'XLim',[1 Truth{1}.K]); set(gca, 'YLim', [-1000,1000]);

% CIS-1
figure(3);
hold on;
subplot(211); box on;

for i = 1:Truth{1}.total_tracks
    hline1 = line(X_track{i}(1,:),X_track{i}(2,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0 * ones(1,3));
end

for k = 1:length(TracksConsensusInfoSelector1)
    colorIndex = mod(k,length(colorBoard));
    if colorIndex == 0
        colorIndex = length(colorBoard);
    end
    if size(TracksConsensusInfoSelector1{k}.X,2) > 0
        tempTime = Time(TracksConsensusInfoSelector1{k}.t);
        hline2 = line(tempTime, TracksConsensusInfoSelector1{k}.X(1,:),'LineStyle','-','Marker','.','Markersize',8,'Color',colorBoard{colorIndex});
    end
end

subplot(212);box on;


for i = 1:Truth{1}.total_tracks
    yhline1 = line(X_track{i}(1,:), X_track{i}(4,:), 'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3));
end

for k = 1:length(TracksConsensusInfoSelector1)
    colorIndex = mod(k,length(colorBoard));
    if colorIndex == 0
        colorIndex = length(colorBoard);
    end
    if size(TracksConsensusInfoSelector1{k}.X,2) > 0
        tempTime = Time(TracksConsensusInfoSelector1{k}.t);
        hline2 = line(tempTime, TracksConsensusInfoSelector1{k}.X(3,:),'LineStyle','-','Marker','.','Markersize',8,'Color',colorBoard{colorIndex});
    end
end

subplot(211); xlabel('Time:s'); ylabel('x:m');
set(gca, 'XLim', [1 Truth{1}.K]); set(gca,'YLim',[-1000,1000]);
legend([hline2 hline1],'Estimation','Ground truth');

subplot(212); xlabel('Time:s'); ylabel('y:m');
set(gca, 'XLim',[1 Truth{1}.K]); set(gca, 'YLim', [-1000,1000]);

% CIS-2
figure(4);
hold on;
subplot(211); box on;

for i = 1:Truth{1}.total_tracks
    hline1 = line(X_track{i}(1,:),X_track{i}(2,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0 * ones(1,3));
end

for k = 1:length(TracksConsensusInfoSelector2)
    colorIndex = mod(k,length(colorBoard));
    if colorIndex == 0
        colorIndex = length(colorBoard);
    end
    if size(TracksConsensusInfoSelector2{k}.X,2) > 0
        tempTime = Time(TracksConsensusInfoSelector2{k}.t);
        hline2 = line(tempTime, TracksConsensusInfoSelector2{k}.X(1,:),'LineStyle','-','Marker','.','Markersize',8,'Color',colorBoard{colorIndex});
    end
end

subplot(212);box on;


for i = 1:Truth{1}.total_tracks
    yhline1 = line(X_track{i}(1,:), X_track{i}(4,:), 'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3));
end

for k = 1:length(TracksConsensusInfoSelector2)
    colorIndex = mod(k,length(colorBoard));
    if colorIndex == 0
        colorIndex = length(colorBoard);
    end
    if size(TracksConsensusInfoSelector2{k}.X,2) > 0
        tempTime = Time(TracksConsensusInfoSelector2{k}.t);
        hline2 = line(tempTime, TracksConsensusInfoSelector2{k}.X(3,:),'LineStyle','-','Marker','.','Markersize',8,'Color',colorBoard{colorIndex});
    end
end

subplot(211); xlabel('Time:s'); ylabel('x:m');
set(gca, 'XLim', [1 Truth{1}.K]); set(gca,'YLim',[-1000,1000]);
legend([hline2 hline1],'Estimation','Ground truth');

subplot(212); xlabel('Time:s'); ylabel('y:m');
set(gca, 'XLim',[1 Truth{1}.K]); set(gca, 'YLim', [-1000,1000]);

% CMIL
figure(5);
hold on;
subplot(211); box on;

for i = 1:Truth{1}.total_tracks
    hline1 = line(X_track{i}(1,:),X_track{i}(2,:),'LineStyle','-','Marker','none','LineWidth',1,'Color',0 * ones(1,3));
end

for k = 1:length(TracksCMIL)
    colorIndex = mod(k,length(colorBoard));
    if colorIndex == 0
        colorIndex = length(colorBoard);
    end
    if size(TracksCMIL{k}.X,2) > 0
        tempTime = Time(TracksCMIL{k}.t);
        hline2 = line(tempTime, TracksCMIL{k}.X(1,:),'LineStyle','-','Marker','.','Markersize',8,'Color',colorBoard{colorIndex});
    end
end

subplot(212);box on;


for i = 1:Truth{1}.total_tracks
    yhline1 = line(X_track{i}(1,:), X_track{i}(4,:), 'LineStyle','-','Marker','none','LineWidth',1,'Color',0*ones(1,3));
end

for k = 1:length(TracksCMIL)
    colorIndex = mod(k,length(colorBoard));
    if colorIndex == 0
        colorIndex = length(colorBoard);
    end
    if size(TracksCMIL{k}.X,2) > 0
        tempTime = Time(TracksCMIL{k}.t);
        hline2 = line(tempTime, TracksCMIL{k}.X(3,:),'LineStyle','-','Marker','.','Markersize',8,'Color',colorBoard{colorIndex});
    end
end

subplot(211); xlabel('Time:s'); ylabel('x:m');
set(gca, 'XLim', [1 Truth{1}.K]); set(gca,'YLim',[-1000,1000]);
legend([hline2 hline1],'Estimation','Ground truth');

subplot(212); xlabel('Time:s'); ylabel('y:m');
set(gca, 'XLim',[1 Truth{1}.K]); set(gca, 'YLim', [-1000,1000]);

