%% s_preprocessMEGRetinotopyData

% This script is in progress but will hopefully replace Barrie's old scripts
% for the main analysis to preprocess the MEG dataset.

% This script relies on the following toolboxes:
% * meg_utils
% * Fieldtrip
% * megdenoise

% All toolboxes can be added with the toolboxtoolbox command:
% tbUse('retmeg') 

%% 0. Define parameters and paths

% Define subject and data path
subject_nr    = 'wlsubj030';
fnameCombined = 'R0942_RetCombined'; % if session had multiple sqd files, we can combine them with another function and load that file
fnameSingle   =  '*_Ret_*';          % .. or just load the individual files
dataPth       = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/';

% Derive other file paths
rawSqdPath = fullfile(dataPth, subject_nr, 'raw');
paramFilePth = fullfile(rawSqdPath, 'R0942_MegRet_9.11.18', 'behavior');
stimFilePth = fullfile(rawSqdPath, 'R0942_MegRet_9.11.18', 'stimFiles');

% Make 'processed' folder to save time series
savePth = fullfile(dataPth, subject_nr, 'processed');
if ~exist(savePth, 'dir'); mkdir(savePth); end

 % MEG Channel information:
triggerChan    = 161:168;
photoDiodeChan = 192;
dataChan       = 1:157;
fs             = 1000; % Sample rate (Hz)

% Preprocessing options:
verbose       = true;
doFiltering   = true;
doDenoise     = true;
doSaveData    = true;

% Go to dataPth
curDir = pwd;
cd(dataPth);

%% 1. Load MEG data

if verbose; sprintf('(%s) Load sqd data...\n', mfilename); end
ts = meg_load_sqd_data(rawSqdPath, fnameSingle);
       

%% 2. Get Triggers
    
if verbose; sprintf('(%s) Get triggers from data...\n', mfilename); end

% Trigger number legend:
% 1-8 = barsweep orientations
% 20  = blank
% 10  = blink
% 

triggers.ts = meg_fix_triggers(ts(triggerChan,:)'); % (ts should be time x chan, function from meg_utils)
triggers.timing = find(triggers.ts);

% remove first trigger (not sure why this one is here)
triggers.ts(triggers.timing(1)) = 0;
triggers.timing = find(triggers.ts);

medianTriggerLength = median(diff(triggers.timing));
if verbose
    fprintf('Median trigger length is %d samples\n', medianTriggerLength);
    triggerConditions = unique(triggers.ts);
    for ii = 1:length(triggerConditions(triggerConditions>0))
        fprintf('Condition %02d has \t%d triggers\n', triggerConditions(ii), sum(triggers.ts==triggerConditions(ii)));
    end
        
    % Plot triggers   
    figure; plot(triggers.ts); title('Triggers'); xlabel('Time (ms)'); ylabel('trigger number');
end


%% 3. Epoching

% Epoch information about all epochs
triggers.stimConditions       = triggers.ts(triggers.ts>0);
totalEpochs                   = length(triggers.stimConditions);

% Epoch information about stimulus (bar sweep) epochs
triggers.onlyBarStim          = find((triggers.ts>0) & (triggers.ts<10));
totalStimEpochs               = length(triggers.onlyBarStim);

% Derived stimulus information
numOfOrientations             = length(triggerConditions((triggerConditions>0)&(triggerConditions<10))); % note, only 5 orientations for subj040 and 004, thus 140 epochs per run
numRuns                       = length(dir(fullfile(stimFilePth,'*stimulus*')));
numOfEpochsPerOrientation     = totalStimEpochs / numOfOrientations / numRuns;

epochStartEnd                 = [.150 (.150+1.100)]; % s (First 150 ms are blank, one epoch length = 1.100 s)
flickerFreq                   = 10; % Hz

if verbose; sprintf('(%s) Epoch data...\n', mfilename); end
data = getEpochs(ts(dataChan,:), triggers, epochStartEnd, flickerFreq, fs); % time x epochs x channels

% Set blink and blank epochs to NaNs
data(:, triggers.stimConditions==10,:) = NaN;
% data(:, triggers.stimConditions==20,:) = NaN;

if verbose
    t = (1:size(data,1))./fs;
    ft = (0:length(t)-1)/max(t);
    figure; subplot(211);
    plot(t,nanmean(data(:,:,1),2));
    xlabel('Time (s)'); ylabel('Magnetic Flux (Tesla)');
    subplot(212);
    amps = abs(fft(data(:,:,1)))/size(data,1)*2;
    plot(ft,nanmean(amps,2));
    xlim([0 100]); xlabel('Frequency (Hz)'); ylabel('Amplitudes');
end

%% 4. Filter MEG data
if doFiltering

    fParams.fStop = 0.1;    % LowPass filter (Hz)
    fParams.fPass = 1;      % BandPass filter (Hz)
    fParams.aStop = 60;     % Db
    fParams.aPass = 3;      % Passband Ripple (Db)
    fParams.fs    = fs;
    
    if verbose; sprintf('(%s) High pass filter data...\n', mfilename); end
    data = highPassFilterData(data, fParams, verbose);

end
%% 5. Label bad channels / epochs

if verbose; sprintf('(%s) Check for bad channels or epochs in data...\n', mfilename); end

varThreshold        = [0.05 20];
badChannelThreshold = 0.2;
badEpochThreshold   = 0.2;
verbose             = true;

[data, badChannels, badEpochs]  = nppPreprocessData(data, ...
    varThreshold, badChannelThreshold, badEpochThreshold, verbose);

%% 6. Denoise time series

if doDenoise

    if verbose; sprintf('(%s) Denoise data...\n', mfilename); end
    
    % Define denoising options
    opt.verbose             = verbose;
    opt.use3Channels        = false;
    opt.removeFirstEpoch    = false;
    opt.removeMsEpochs      = false;
    opt.pcchoose            = 1.05;  % initial threshold
    opt.npcs2try            = 10;    % max nr of PCs = 10

    % Get 10 Hz evoked signal (fft power)
    evokedfun        = @(x)mprfDenoiseEvalFun(x,[flickerFreq, flickerFreq] ,fs);
    designConditions = triggers.stimConditions;
    designConditions(designConditions > 9)=0;
    designMatrix     = conditions2design(designConditions);

    dataToDenoise = permute(data, [3,1,2]);
    [results,evalout,denoised_spec, denoised_data] = ...
        denoisedata(designMatrix(~badEpochs,:),dataToDenoise(~badChannels, :, ~badEpochs), evokedfun, evokedfun, opt);

    if verbose
        figure; plot(results.origmodel.r2,results.finalmodel.r2,'k.'); hold on;
        plot([0 max([results.origmodel.r2,results.finalmodel.r2])],[0 max([results.origmodel.r2,results.finalmodel.r2])],'r-');
        ylabel('Final R2');
        xlabel('Orig R2');
        axis square;
        

        amps = abs(fft(denoised_data{1}(1,:,:),[],2))/size(data,2)*2;
        figure; plot(ft,nanmean(squeeze(amps),2));
        xlim([0 100]); xlabel('Frequency (Hz)'); ylabel('Amplitudes');
    end

    dataDenoised = NaN(size(denoised_data{1}));
    dataDenoised(~badChannels, :, ~badEpochs) = denoised_data{1};
end

%% 7. Reblock and save prepreocessed timeseries

% we split the number of runs in half to match with previous datasets

if doSaveData

    if verbose; sprintf('(%s) Save data...\n', mfilename); end
    
    % to do reshape back to time x epochs x channels x runs (maybe split up
    % the runs to get same as old data sets)
    clear dataBlocked;
    numTrigPerBlock = (numOfEpochsPerOrientation + (3*(numOfEpochsPerOrientation+5)));
    startEpoch = 1;
    figure; plot(triggers.stimConditions); hold all;
    for n = 1:numRuns*2
        
        if n == 1
            theseEpochs = startEpoch:(startEpoch+numTrigPerBlock-1);
        elseif (mod(n,2)==1)
            startEpoch = theseEpochs(end)+1;
            theseEpochs = startEpoch:(startEpoch+numTrigPerBlock-1);   
        else
            startEpoch = theseEpochs(end)+6;
            theseEpochs = startEpoch:(startEpoch+numTrigPerBlock-1);
        end
        
        dataBlocked(:,:,:,n) = dataDenoised(:,:, theseEpochs);
        plot(theseEpochs,triggers.stimConditions(theseEpochs),'r:', 'LineWidth', 4);
%         waitforbuttonpress;
    end   
    
    data = permute(dataBlocked, [2, 3, 4, 1]); % channels x time x epochs x blocks --> time x epochs x blocks x channels 
    save(fullfile(savePth, 'MEG_timeseries.mat'), 'data', '-v7.3')

end






