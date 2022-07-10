%Analyze data from recordings of V1 cells responding to oriented
%gratings presented at various angles
clear all
close all
%LOAD IN DATA
load Sur_Orientation_SpikeData
spikes = double(spikes);
%SET UP USEFUL VARIABLES
NumControls = 2; %number of control experiments with no grating
dt = 1; %spacing between sampled time points [ms]
t_On = 0; %time stimulus turns on [ms]
t_Move = 500; %time stimulus begins moving [ms]
t_Off = 2500; %time stimulus turns off [ms]
NumAngles = size(spikes,1) - NumControls %number of angles tested, equally spaced;
%last 2 sets of recordings are controls
NumTimePoints = size(spikes,2) %number of time points; time was sampled every 1 ms
NumTrials = size(spikes,3) %number of trials performed at each angle
t_vect = t_On:dt:(NumTimePoints-1)*dt; %time vector for each trial
ThisOrientation = 13; %element index of orientation we are currently analyzing
%PLOT RASTERS FOR ONE PARTICULAR ANGLE
figure(1)
for trial=1:NumTrials
plot(t_vect,trial*spikes(ThisOrientation,:,trial),'+')
hold on
end
xlabel('time (ms)')
ylabel('trial number')
axis([0 3500 0.5 NumTrials])
hold off
%PLOT AVERAGED DATA (PSTH) FOR ONE PARTICULAR ANGLE
TrialSum_vect = sum(spikes(ThisOrientation,:,:),3); %adds up spike trains across trials
SmoothingWidth = 61; %smooths data over +/- [(this # - 1)/2] bins
TrialSum_smooth_vect = smooth(TrialSum_vect,SmoothingWidth);
figure(2)
plot(t_vect,1000*TrialSum_smooth_vect/(NumTrials*dt))
xlabel('time (ms)')
ylabel('Average firing rate (Hz)')


%COMPUTE FANO FACTOR
TStartCount = 600; %time to start computing average
TEndCount = 2500; %time to end computing average
%next line gives number of counts for each trial
NumCounts_vect = sum(spikes(ThisOrientation,((TStartCount+1)/dt):(TEndCount/dt),:),2);
FanoFactor = (std(NumCounts_vect)^2)/mean(NumCounts_vect)
%COMPUTE CV_isi
SpikeTimes_vect = dt*find(abs(spikes(ThisOrientation,TStartCount:TEndCount,1)-1) < 0.00000001);
isi_vect = diff(SpikeTimes_vect);
%plot ISI histogram
figure(3)
hist(isi_vect,8)
xlabel('isi (ms)')
ylabel('Number of occurrences')


% %compute CV_isi
% mean_isi = mean(isi_vect)
% std_isi = std(isi_vect)
% CV_isi = mean_isi/std_isi
% %COMPUTE AVE RATE FOR 1 TRIAL
% TotalCounts = sum(NumCounts_vect) %total spikes across time and trials
% AveRate = 1000*TotalCounts/((TEndCount-TStartCount)*NumTrials) %average across time and trials
% %COMPUTE TUNING CURVES
% %next line gives matrix of counts across trials for various orientations and times
% TuningCurveCounts_matrix = sum(spikes(:,((TStartCount+1)/dt):(TEndCount/dt),:),3);
% TuningCurveCounts_vect = sum(TuningCurveCounts_matrix,2) %sum over time dimension
% TuningCurve_AveRate_vect = 1000*TuningCurveCounts_vect/((TEndCount-TStartCount)*NumTrials)
% %plot tuning curve
% DeltaTheta = 360/NumAngles
% Orientation_vect = 0:DeltaTheta:(NumAngles-1)*DeltaTheta;
% figure(4)
% plot(Orientation_vect,TuningCurve_AveRate_vect(1:NumAngles))
% xlabel('orientation (deg)')
% ylabel('Ave firing rate (Hz)')
% Control_AveRate = mean(TuningCurve_AveRate_vect(17:18))
% %plot difference from control average rate
% TuningCurve_RateDiff_vect = TuningCurve_AveRate_vect - Control_AveRate;
% figure(5)
% plot(Orientation_vect,TuningCurve_RateDiff_vect(1:NumAngles))
% xlabel('orientation (deg)')
% ylabel('Ave firing rate relative to control(Hz)')
% grid on