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
subject       = 'wlsubj058';
fnameSingle   =  '*Ret*';          % case sensitive!
dataPth       = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/';

% Derive other file paths
rawSqdPath = fullfile(dataPth, subject, 'raw');
paramFilePth = fullfile(dataPth, subject,  'paramFiles');
stimFilePth = fullfile(dataPth, subject, 'stimFiles');

% Make 'processed' folder to save time series
savePth = fullfile(dataPth, subject, 'processed');
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
verbose       = true;

% Go to dataPth
curDir = pwd;
cd(dataPth);

%% 1. Load MEG data

if verbose; sprintf('(%s) Load sqd data...\n', mfilename); end
[ts, meg_files] = meg_load_sqd_data(rawSqdPath, fnameSingle);
       

%% 2. Get Triggers
    
if verbose; sprintf('(%s) Get triggers from data...\n', mfilename); end

% Trigger number legend:
% 1-8 = barsweep orientations
% 20  = blank
% 10  = blink
% 

if strcmp('wlsubj058',subject) % subject got some trigger missings during experiment, thus needs its own function
    triggers.ts = meg_fix_triggers_wlsubj058(ts, triggerChan);
elseif strcmp('wl_subj040', subject)  % subject got random extra triggers and needs its own function
    triggers.ts = meg_fix_triggers_wlsubj040(ts(triggerChan,:)');
else
    triggers.ts = meg_fix_triggers(ts(triggerChan,:)'); % (ts should be time x chan, function from meg_utils)
end

triggers.timing = find(triggers.ts);

% For subject wlsubj030: remove first trigger (not sure why this one is here)
if strcmp(subject, 'wlsubj030'); triggers.ts(triggers.timing(1)) = 0; triggers.timing = find(triggers.ts); end;
% For subject wlsubj004: remove half run
if strcmp(subject, 'wl_subj004'); triggers.ts(697111:729608) = 0; triggers.timing = find(triggers.ts); end


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

% save stimulus conditions
save(fullfile(savePth, 'megStimConditions.mat'), 'triggers')


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
[data, startOfRun] = getEpochs(ts(dataChan,:), triggers, epochStartEnd, flickerFreq, fs); % time x epochs x channels

% get size of data before removing epochs
sz = size(data);

% Set blink epochs to NaNs, leave blanks as is
data(:, triggers.stimConditions==10,:) = NaN;   % 10: Blink blocks
% data(:, triggers.stimConditions==20,:) = NaN; % 20: Blank blocks

% Set the first epoch of the first bar sweep to NaN ([EK]: not sure why we
% do this? Maybe consider not doing this since we already remove the first
% 150 ms of every epoch?)
data(:, startOfRun, :) = NaN;

% Plot a single channel to check data
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
    fParams.aStop = 60;     % ?? Amplitude of lowpass filter(Db)
    fParams.aPass = 3;      % ?? Amplitude band pass Ripple (Db)
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
    if strcmp(subject, 'wlsubj068') || strcmp(subject, 'wlsubj058') || strcmp(subject, 'wl_subj004') || strcmp(subject, 'wl_subj040')
        designConditions(designConditions==3) = 2;
        designConditions(designConditions==4) = 3;
        designConditions(designConditions==6) = 4;
        designConditions(designConditions==7) = 5;
    end
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

    dataDenoised = NaN(sz(3),sz(1),sz(2));
    dataDenoised(~badChannels, :, ~badEpochs) = denoised_data{1};
end

%% 7. Reblock and save prepreocessed timeseries

% we split the number of runs in half to match with previous datasets

if doSaveData
    numRuns = 19;
    % to do reshape back to time x epochs x channels x runs
    if verbose; sprintf('(%s) Save data...\n', mfilename); end
     
    % SPLIT UP IN blocks   (i.e. number of RUNS)
    dataBlocked =  NaN(sz(3), sz(1), sz(2)/numRuns, numRuns);
    
     % plot all triggers
    figure; plot(triggers.stimConditions, 'LineWidth', 2); xlabel('Time (timepoints)'); ylabel('Trigger num'); hold all;
    startOfRun = 1:140:(19*140);
    for n = 1:numRuns
       
        % Get first epoch
        startEpoch = startOfRun(n);

        % Get last epoch depending on the run number
        if n < numRuns
            lastEpoch  = startOfRun(n+1)-1;
            theseEpochs = startEpoch:lastEpoch;
        else
            lastEpoch  = size(dataBlocked,3);
            theseEpochs = startEpoch:(startEpoch+lastEpoch-1);
        end
        
        % Get data and put in dataBlocked variable
        dataBlocked(:,:,:,n) = dataDenoised(:,:, theseEpochs);

        % Mark those that are blocked
        plot(theseEpochs,triggers.stimConditions(theseEpochs),'r:', 'LineWidth', 4);
    end
    
    clear data;
    data.data = permute(dataBlocked, [2, 3, 4, 1]); % channels x time x epochs x blocks --> time x epochs x blocks x channels
    save(fullfile(savePth, 'epoched_data_hp_preproc_denoised.mat'), 'data', '-v7.3')
end






