%% makeFigureXX_eyeMovements_xyPosition
% Script to plot horizontal and vertical eye movement for MEG scan run for
% MEG retinotopy project. 
% 
% need scripts from noisepoolPCADenoise project - 
% addpath(genpath('~/noisepoolPCADenoise/external/eyetracking'));
% addpath(genpath('/Volumes/server/Projects/MEG/Eyetracking_scripts/'));

subjID   = 'wlsubj081'; % For MEG data choose from wlsubj039, wlsubj058, wlsubj068, wlsubj070, wlsubj081, wlsubj106, wlsubj109, wlsubj111
modality = 'MEG'; % for now, we can only process data from 'MEG' and not from 'MRI' (Note: not all subjects have both datasets)

% Load paths with data files for this subject
dirPth = loadPaths(subjID);

% Go back to root
cd(mprf_rootPath)

% Set options
opt = getOpts('saveFig',1,'verbose',1); % see getOpts function for more options

% Load eye link mat file
load(fullfile(dirPth.(lower(modality)).eyePth,sprintf('%s_eyelink.mat',modality)),'eyd');

% Load stimulus triggers 
load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'megStimConditions.mat'), 'triggers');

% Get start time in eyelink data
for ii = 1:length(eyd.messages)
    if strfind(eyd.messages(ii).message, 'MEG Trigger')
        startTime = ii; break;
    end
end

% Get xy data and messages
endTime     = 0; % number of messages to omit from end (if any)
timeLims    = [eyd.messages(startTime).time(1), eyd.messages(end-endTime).time];
blinkWinSec = [0.2 0.35];
s           = msGetEyeData(eyd,timeLims,blinkWinSec);

% Find start of run (20 is blink block)
blinkTrigNr   = 20;
blinkIndex    = find(triggers.stimConditions==blinkTrigNr);
stimIndex     = find(triggers.stimConditions~=blinkTrigNr);

runStartIndex    = blinkIndex(1:2:end);
diffNrOfEpochs   = diff(runStartIndex);
nrOfEpochsPerRun = diffNrOfEpochs(1); % note, after the first blank epoch


% Define other run details
epochDur          = 1300; % 1300 ms with a sampling rate of 1000 Hz
stimRad           = 10;   % deg
nBarPassesPerRun  = 5;
runs              = 19;

% Find the messages with run start
blinkMessageIndex      = find(strcmp({eyd.messages.message}, sprintf('MEG Trigger: %d',blinkTrigNr)) == 1 ); % every run starts with message - MEG trigger: 20
runStartMessageIndex   = blinkMessageIndex(1:2:end);

allTimePoints = [];
for runIdx = 1:length(runStartMessageIndex)
    
   if mod(runIdx,6) % after 5 bar passes, 6th run there is eye movement (there was a break ?) - so that was excluded
        
       startOfRun = runStartMessageIndex(runIdx);
       
       startOfRun = startOfRun+2; % skip the two blank blocks
       
       epochsToInclude  = startOfRun:(startOfRun+nrOfEpochsPerRun-3); % 27 epochs (excluding 2 blinks, but including 3 blanks and 22 bar pass)
        
       epochsToInclude_time = [eyd.messages(epochsToInclude).time] - s.timeRaw(1); % time is calculated with respect to a reference time 
        
       thisRunTimePoints = epochsToInclude_time'+(0:epochDur-1); % there are 1300 points between consecutive epochsToInclude_time values because every epoch was 1300 ms long

       allTimePoints = [allTimePoints; thisRunTimePoints(:)]; % total should be 1300*25*5*19 = 3087500
    end
    
end

xPos = s.xyPos(allTimePoints,1);
yPos = s.xyPos(allTimePoints,2);

% Get mean, median and covariance matrix
traceMean   = nanmean([xPos, yPos]);
traceMedian = nanmedian([xPos, yPos]);
traceCov    = cov([xPos(~isnan(xPos)), yPos(~isnan(yPos))]);

% plot horizontal vs vertical eyemovement positions
fH = figure(1); clf; 
set(gcf,'Color', 'w', 'Position',[31, 252, 1067, 996]);  hold on;
% Plot eye traces
plot(xPos,yPos, '.','MarkerSize',1, 'Color', [.7 .7 .7]);
% Plot 95% confidence ellipse
ax = error_ellipse(traceCov,traceMean,'conf',0.95,'color','k');
% Plot median
plot(traceMean(1),traceMean(2), 'r+', 'markersize',5);

% Scale axes
xlim([-stimRad stimRad]);
ylim([-stimRad stimRad]);
xlabel('X Position (deg of visual angle)');
ylabel('Y Position (deg of visual angle)');
grid on; box off;
axis square; 
title(sprintf('Subject %s - Fixation', subjID))
set(gca,'FontUnits','centimeters', 'FontSize', 0.5,'TickDir','out','LineWidth',2);


% save the figure 
if opt.saveFig
    [pth, ~] = fileparts(dirPth.model.saveFigPth);
    saveDir = fullfile(pth, 'finalfig', 'figureXX_EyeFixation');
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    fprintf('\n(%s): Saving eye movements figure  in %s\n',mfilename, saveDir);

    print(fH, fullfile(saveDir, sprintf('figXX_eyeMovements_%s_%s', subjID,lower(modality))), '-dpng');
    print(fH, fullfile(saveDir, sprintf('figXX_eyeMovements_%s_%s', subjID,lower(modality))), '-depsc');

end






