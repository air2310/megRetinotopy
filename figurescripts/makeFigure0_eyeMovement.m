% makeFigure0_eyeMovement.m
% Script to plot horizontal and vertical eye movement for MEG scan run for
% MEG retinotopy project. 
% 
% need scripts from noisepoolPCADenoise project - 
% addpath(genpath('~/noisepoolPCADenoise/external/eyetracking'));

subjID = 'wlsubj081';
% Load paths with data files for this subject
dirPth = loadPaths(subjID);

% Go back to root
cd(mprf_rootPath)

% Set options
opt = getOpts('saveFig',1,'verbose',1); % see getOpts function for more options

eyeFilePath = fullfile(dirPth.meg.eyePth); 

eyeFileName = 'MEG_eyelink.mat';
eyeFile     = fullfile(eyeFilePath,eyeFileName);

% load eye data 
load(eyeFile,'eyd');
s = msGetEyeData(eyd);

% stimulus parameters 
blinkEpochs  = 2;
blankEpochs  = 3;
barEpochs    = 22;
runEpochs    = (blinkEpochs + blankEpochs + barEpochs);
epochDur     = 1300; % 1300 ms with a sampling rate of 1000 Hz
stimRad      = 10; % 10 degrees
barPass      = 5;
runs         = 19;

% find the epoch start
blinkIndex      = find(strcmp({eyd.messages.message}, 'MEG Trigger: 20') == 1 ); % every run starts with message - MEG trigger: 20
runStartIndex   = blinkIndex(1:2:end);

Epoch_time_all = [];
%figure
for runIdx = 1:size(runStartIndex,2)
    
    if mod(runIdx,6) % after 5 bar passes, 6th run there is eye movement (there was a break ?) - so that was excluded
        
        StartOfEpoch = (runStartIndex(runIdx):runStartIndex(runIdx)+runEpochs); % 28 epochs including blink, blanks and bar pass
        
        StartOfEpoch_time = [eyd.messages(StartOfEpoch).time] - s.timeRaw(1); % time is calculated with respect to a reference time 
        
        Epoch_time = StartOfEpoch_time(1):StartOfEpoch_time(end)-1; % there are 1300 points between consecutive StartOfEpoch_time values because every epoch was 1300 ms long
        
        Epoch_time_all = [Epoch_time_all Epoch_time]; % combine all the time points from 19 runs (each with 5 bar passes) - 6 (5 + 1(break after a run)) * 19 = 114 epochs  
        
    end
    
end


% plot horizontal vs vertical eyemovement positions
fH = figure(1); clf;
set(gcf,'position',[66,1,1855,1001]);
plot(s.xyPos(Epoch_time_all,1),s.xyPos(Epoch_time_all,2),'color',[0 0 0]);
xlim([-stimRad stimRad]);
ylim([-stimRad stimRad]);
xlabel('degrees of visual angle');
ylabel('degrees of visual angle');
grid on;
axis square;
set(gca,'FontUnits','centimeters', 'FontSize', 0.5,'TickDir','out','LineWidth',3);


% save the figure 
if opt.saveFig
    [pth, ~] = fileparts(dirPth.model.saveFigPth);
    saveDir = fullfile(pth, 'finalfig', 'figure0');
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    fprintf('\n(%s): Saving eye movements figure  in %s\n',mfilename, saveDir);

    print(fH, fullfile(saveDir, sprintf('eyeMovements%s_%s', dirPth.subjID,'meg')), '-dpng');
end






