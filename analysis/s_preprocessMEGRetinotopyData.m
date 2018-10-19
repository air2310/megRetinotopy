%% s_preprocessMEGRetinotopyData

% This script is the main analysis to preprocess the MEG dataset.

% This script relies on the following toolboxes:
% * meg_utils
% * Fieldtrip
% * megdenoise

% all toolboxes can be added with the toolboxtoolbox command tbUse('retmeg') 

%% 0. Define parameters and paths

% Define subject and data path
subject_nr    = 'wlsubj030';
fnameCombined = 'R0942_RetCombined'; % if session had multiple sqd files, we can combine them
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
if verbose; figure; plot(triggers.ts); title('Triggers'); xlabel('Time (ms)'); ylabel('trigger number'); end


%% 3. Epoching

% Epoch information
totalNumEpoch                 = length(triggers.timing);
numStimConditionsWithinOneRun = 211; % Note: 140 for subject004 and 040
numRuns = totalNumEpoch/numStimConditionsWithinOneRun;
epochStartEnd                 = [150 (150+1100)]; % ms (First 150 ms are blank, total amount of timepoints = 1100)
flickerFreq                   = 10; % Hz

if verbose; sprintf('(%s) Epoch data...\n', mfilename); end
data = getEpochs(ts(dataChan,:), triggers, epochStartEnd, flickerFreq, fs);


%% 4. Filter MEG data
if doFiltering

    if verbose; sprintf('(%s) High pass filter data...\n', mfilename); end
    data = mprfFilterDataHighPass(data);

end
%% 5. Label bad channels / epochs

if verbose; sprintf('(%s) Check for bad channels or epochs in data...\n', mfilename); end
[data, badChannels, badEpochs] = mprfPreprocessWrapper(data);


%% 6. Denoise time series

if doDenoise

    if verbose; sprintf('(%s) Denoise data...\n', mfilename); end
    data = mprfDenoiseWrapper(data, stimFilePth);

end
%% 7. Save prepreocessed timeseries
if doSaveData

    if verbose; sprintf('(%s) Save data...\n', mfilename); end
    save(fullfile(savePth, 'MEG_timeseries.mat'), 'epochedDataFilteredCleanDenoised')

end






