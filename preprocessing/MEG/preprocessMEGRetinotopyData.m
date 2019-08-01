function [data, triggers, opt] = preprocessMEGRetinotopyData(subjID, dirPth, opt)
% This function contains main analysis to preprocess subject's MEG data,
% from the MEG Retinotopy project.
%
%       data = preprocessMEGRetinotopyData(subjID, dirPth, opt)
%
% INPUTS
%   subjID :      subject nr (string)
%   dirPth :      paths locating subject's data and files (struct, see loadPaths.m)
%   opt    :      options for preprocessing (struct)
%                     - verbose     :  print text in command window and
%                                      show debug figures
%                     - doFiltering :  low pass filter
%                     - doDenoise   :  GLM denoise for MEG
%                     - doSaveData  :  Save data and conditions as matfiles
%                     - saveFig     :  Save debug figures, only if verbose = true
%                     - removeStartOfRunEpoch :  remove the first epoch of
%                                                the bar sweep
% OUTPUTS:
%   preprocData    :     preprocessed data (time x epochs x channels x runs)
%
% WORKFLOW:
%   0. Define parameters and paths
%   1. Load MEG data
%   2. Get Triggers
%   3. Epoching
%   4. Filter MEG data
%   5. Label bad channels / epochs
%   6. Denoise time series
%   7. Reblock and save prepreocessed timeseries
%
% DEPENDENCIES:
% This script relies on the following toolboxes:
%   * meg_utils
%   * Fieldtrip
%   * megdenoise
%
% All dependent toolboxes can be added with the ToolboxToolbox command:
%   tbUse('retmeg')
%
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019

%% 0. Define parameters and paths

% Make folder to save figure
if (opt.saveFig && ~exist(dirPth.meg.saveFigPth, 'dir'))
    mkdir(dirPth.meg.saveFigPth);
end

% Define and create 'processed' folder to save time series
if opt.removeStartOfRunEpoch
    saveDataPth = fullfile(dirPth.meg.processedDataPth, 'firstEpochRemoved');
else
    saveDataPth = fullfile(dirPth.meg.processedDataPth, 'allEpochs');
end

if ~exist(saveDataPth, 'dir'); mkdir(saveDataPth); end

% Go to dataPth
curDir = pwd;
cd(fullfile(dirPth.meg.dataPth,subjID));

%% 1. Load MEG data

if opt.verbose; sprintf('(%s) Load sqd data...\n', mfilename); end
[ts, megFiles] = meg_load_sqd_data(dirPth.meg.rawSqdPth, '*Ret*'); %#ok<ASGLU>


%% 2. Get Triggers

% Trigger number legend:
% 20  = blink - every run starts with 2 blink TRs
% 10  = blank - followed by 3 blank TRs,
% followed by stimulus:
% 1-8 = barsweep orientations, there are only 5 orientations in the MEG
% experiment with the following trigger nr's: 1, 3, 4, 6, 7.
%
% Some subjects have triggers missing (wlsubj058), or extra random triggers
% (wlsubj040) which cannot be filtered out by our general meg_fix_triggers
% function. Therefore, these subjects have their own (modified) version of
% that function. Output of triggers.ts is time x chan (function from meg_utils)

if opt.verbose; fprintf('(%s) Get triggers from data...\n', mfilename); end

% Get trigger time series (same length as MEG ts) and compute triggertiming
switch subjID
    case 'wlsubj004'
        triggers.ts = meg_fix_triggers(ts(opt.triggerChan,:)');
        triggers.ts(697111:729608) = 0;      % wlsubj004: remove half run (to do: make it a general statement, not hardcoded)
        triggers.timing = find(triggers.ts);
    case 'wlsubj030'
        triggers.ts = meg_fix_triggers(ts(opt.triggerChan,:)');
        triggers.timing = find(triggers.ts);
        triggers.ts(triggers.timing(1)) = 0; % wlsubj030: remove first trigger (not sure why this one is here)
        triggers.timing = find(triggers.ts);
    case 'wlsubj040'
        triggers.ts = meg_fix_triggers_wlsubj040(ts(opt.triggerChan,:)');
        triggers.timing = find(triggers.ts);
    case 'wlsubj058'
        triggers.ts = meg_fix_triggers_wlsubj058(ts,opt.triggerChan);
        triggers.timing = find(triggers.ts);
    case 'wlsubj068'
        triggers.ts = meg_fix_triggers(ts(opt.triggerChan,:)');
        triggers.timing = find(triggers.ts);
    case 'wlsubj039'
        triggers.ts = meg_fix_triggers(ts(opt.triggerChan,:)');
        triggers.timing = find(triggers.ts);
    otherwise
        error('(%s) Can''t find subject ID', mfilename)
end

% Check the median trigger length (should be 1300 ms, note that after epoching this will be 1100 ms)
medianTriggerLength = median(diff(triggers.timing));

if opt.verbose
    fprintf('(%s) Median trigger length is %d samples\n', mfilename, medianTriggerLength);
    triggerConditions = unique(triggers.ts);
    for ii = 1:length(triggerConditions(triggerConditions>0))+1
        fprintf('Condition %02d has \t%d triggers\n', triggerConditions(ii), sum(triggers.ts==triggerConditions(ii)));
    end
    
    % DEBUG: Plot triggers
    figure; plot(triggers.ts); title('Triggers'); xlabel('Time (ms)'); ylabel('trigger number');
    set(gca, 'TickDir', 'out', 'FontSize', 14);
    if opt.saveFig
        print(gcf, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_triggers', subjID)))
    end
end


%% 3. Epoching

% Epoch information about all epochs
triggers.stimConditions       = triggers.ts(triggers.ts>0);
totalEpochs                   = length(triggers.stimConditions);

% Save stimulus conditions (if requested)
if opt.doSaveData; save(fullfile(saveDataPth, 'megStimConditions.mat'), 'triggers'); end

% Epoch information about stimulus (bar sweep) epochs
triggers.onlyBarStim          = find((triggers.ts>0) & (triggers.ts<10));
totalStimEpochs               = length(triggers.onlyBarStim);

% Derived stimulus information
numOfOrientations             = length(unique(triggers.stimConditions((triggers.stimConditions>0)&(triggers.stimConditions<10)))); % note, only 5 orientations in MEG
numRuns                       = length(dir(fullfile(dirPth.meg.stimFilePth,'MEG_retinotopy_stimulus_run_*')));
numOfEpochsPerRun             = totalEpochs / numRuns;
numOfEpochsPerOrientation     = totalStimEpochs / numOfOrientations / numRuns;

% Print info (if requested)
if opt.verbose
    fprintf('(%s) Epoch data...\n', mfilename);
    fprintf('(%s) Number of epochs: %d, with a total of %d stimulus epochs\n', mfilename, totalEpochs, totalStimEpochs);
    fprintf('(%s) Number of orientations: %d\n', mfilename, numOfOrientations);
    fprintf('(%s) Number of runs: %d\n', mfilename, numRuns);
    fprintf('(%s) Number of epochs per stimulus orientation: %d\n', mfilename, numOfEpochsPerOrientation);
end
% Epoch matrix dimensions: time x epochs x channels
[data, startOfRun] = getEpochs(ts(opt.dataChan,:), triggers, opt.epochStartEnd, opt.flickerFreq, opt.fs, subjID);

if length(startOfRun) > numRuns+1
    warning('(%s) Variable "startOfRun" has more values (%d) than actual number of runs (%d).. resetting "startOfRuns".. \n', mfilename, length(startOfRun), numRuns)
    startOfRun = 1:numOfEpochsPerRun:totalEpochs;
end

% Get size of data before removing epochs
sz = size(data);

% Set blink epochs to NaNs
data(:, triggers.stimConditions==20,:) = NaN;     % 20: Blink blocks

% If requested: set the first epoch of the first bar sweep to NaN
if opt.removeStartOfRunEpoch
    numBlinkBlanks = 5;
    if strcmp(subjID,'wlsubj030')
        beginOfSweep = 1:(numOfEpochsPerOrientation+numBlinkBlanks):numOfEpochsPerRun;
        toRemove = startOfRun+beginOfSweep-1;
    else
        beginOfSweep = (numBlinkBlanks+1):(numOfEpochsPerOrientation+numBlinkBlanks):numOfEpochsPerRun;
        toRemove = startOfRun+beginOfSweep-1;
    end
    
    if opt.verbose
        figure; plot(triggers.stimConditions);
        hold on; plot(toRemove,triggers.stimConditions(toRemove), 'rx')
    end
    data(:, toRemove, :) = NaN;
end

% DEBUG: Plot a single channel to check data
if opt.verbose
    t = (0:size(data,1)-1)./opt.fs;
    ft = (0:length(t)-1)/max(t);
    freqIdx = mprfFreq2Index(size(data,1), opt.flickerFreq, opt.fs);
    stimBarOnly = triggers.stimConditions<10;
    blankOnly = triggers.stimConditions==10;
    
    figure; set(gcf, 'Position', [1000, 260, 887, 1078])
    for chan = 1:size(data,3)
        cla
        subplot(311);
        plot(t,nanmean(data(:,stimBarOnly,chan),2));
        xlabel('Time (s)'); ylabel('Magnetic Flux (Tesla)');
        title(sprintf('Mean timeseries: Sensor %d', chan));
        set(gca, 'TickDir', 'out', 'FontSize', 14);
        drawnow;
        
        subplot(312);cla;
        ampsStim = abs(fft(data(:,stimBarOnly,chan)))/size(data,1)*2;
        ampsBlank = abs(fft(data(:,blankOnly,chan)))/size(data,1)*2;
        plot(ft,nanmean(ampsStim,2)); hold on;
        plot(ft,nanmean(ampsBlank,2)); legend({'Stim', 'Blank'});
        plot([opt.flickerFreq opt.flickerFreq], [min(nanmean(ampsStim,2)) max(nanmean(ampsStim,2))])
        plot([60 60], [min(nanmean(ampsStim,2)) max(nanmean(ampsStim,2))])
        plot([120 120], [min(nanmean(ampsStim,2)) max(nanmean(ampsStim,2))])
        xlim([1 150]); ylim([10^-15 10^-12]);
        xlabel('Frequency (Hz)'); ylabel('FT Amplitudes (Tesla)');
        title(sprintf('Mean amplitudes incoherent spectrum: Sensor %d', chan))
        set(gca, 'TickDir', 'out', 'FontSize', 14, 'XScale', 'log', 'YScale', 'log');
        drawnow;
        
        subplot(313);cla;
        meanStim = nanmean(data(:,stimBarOnly,chan),2);
        meanBlank = nanmean(data(:,blankOnly,chan),2);
        
        ampsStim = abs(fft(meanStim))/size(data,1)*2;
        ampsBlank = abs(fft(meanBlank))/size(data,1)*2;
        plot(ft,ampsStim); hold on;
        plot(ft,ampsBlank); legend({'Stim', 'Blank'});
        plot([opt.flickerFreq opt.flickerFreq], [0 max(nanmean(ampsStim,2))])
        plot([60 60], [min(ampsStim) max(ampsStim)])
        plot([120 120], [min(ampsStim) max(ampsStim)])
        xlim([1 150]); ylim([10^-16 10^-13]);
        xlabel('Frequency (Hz)'); ylabel('FT Amplitudes (Tesla)');
        title(sprintf('Mean amplitudes coherent spectrum: Sensor %d', chan))
        set(gca, 'TickDir', 'out', 'FontSize', 14, 'XScale', 'linear', 'YScale', 'linear');
        drawnow;
        
        if opt.saveFig
            if ~exist(fullfile(dirPth.meg.saveFigPth,'predenoise_timeseries'), 'dir')
                mkdir(fullfile(dirPth.meg.saveFigPth,'predenoise_timeseries'));
            end
            print(gcf, '-dpng', fullfile(dirPth.meg.saveFigPth,'predenoise_timeseries', ...
                sprintf('%s_singleSensor%d_timeseries', subjID, chan)))
        end
        
    end
end

%% 4. Filter MEG data

if opt.doFiltering
    if opt.verbose; sprintf('(%s) High pass filter data...\n', mfilename); end
    
    fParams.fStop = 0.1;    % LowPass filter (Hz)
    fParams.fPass = 1;      % BandPass filter (Hz)
    fParams.aStop = 60;     % Amplitude of lowpass filter (dB)
    fParams.aPass = 3;      % Amplitude band pass Ripple (dB)
    fParams.fs    = opt.fs;     % Sample rate (Hz)
    
    
    data = highPassFilterData(data, fParams, opt.verbose);
end

%% 5. Label bad channels / epochs

if opt.verbose; fprintf('(%s) Check for bad channels or epochs in data...\n', mfilename); end

[data, badChannels, badEpochs]  = nppPreprocessData(data, ...
    opt.varThreshold, opt.badChannelThreshold, opt.badEpochThreshold, opt.verbose);

if strcmp(subjID, 'wlsubj040')
    badChannels([98, 42]) = 1;
end

opt.badChannels = badChannels;
opt.badEpochs   = badEpochs;

if opt.saveFig; print(gcf, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_badChannelsEpochs', subjID))); end

%% 6. Denoise time series

if opt.doDenoise
    if opt.verbose; sprintf('(%s) Denoise data...\n', mfilename); end
    
    % Define denoising options
    opt.use3Channels        = false;
    opt.removeFirstEpoch    = false;
    opt.removeMsEpochs      = false;
    opt.pcchoose            = 1.05;  % initial threshold
    opt.npcs2try            = 10;    % max nr of PCs = 10
    
    % Get 10 Hz evoked signal (fft power) to define noisepool and result
    switch subjID
        case {'wlsubj039', 'wlsubj040'} % get noisepool from coherent spectrum
            freqIdx = mprfFreq2Index(size(data,1), opt.flickerFreq, opt.fs);
            
            epochsToKeep = triggers.stimConditions<10;
            epochsToKeep = epochsToKeep(~badEpochs);
            meanData = squeeze(nanmean(data(:,epochsToKeep,~badChannels),2));
            cohSpectrum = abs(fft(meanData))/size(meanData,1)*2;
            response10Hz = cohSpectrum(freqIdx,:);
            [val, idx] = sort(response10Hz, 'ascend');
            noisePool = zeros(1,length(idx)); noisePool(idx(1:75))=1;
            if opt.verbose
                figure; megPlotMap(to157chan(noisePool,~badChannels,'nans'),  ...
                    [0 1], [], [],[],[],[],'interpolation', 'nearest'); end
            noisePoolFun = logical(noisePool);
            
            % Define function to get results
            evokedfun = @(x)mprfDenoiseEvalFun(x,[opt.flickerFreq, opt.flickerFreq] ,opt.fs);

        otherwise % evokedfun uses incoherent spectrum SSVEP
            % Define evoked signal function to get results (SSVEF from incoherent spectrum)
            evokedfun = @(x)mprfDenoiseEvalFun(x,[opt.flickerFreq, opt.flickerFreq] ,opt.fs);
            
            % Define function to get noisepool (same in this case as evokedfun)
            noisePoolFun = evokedfun;
    end
    
    % Get design matrix
    designConditions = zeros(size(triggers.stimConditions));
    uniqueStimConds = unique(triggers.stimConditions);
    uniqueStimConds = uniqueStimConds(uniqueStimConds<10);
    
    for ii = 1:length(uniqueStimConds)
        designConditions(triggers.stimConditions==uniqueStimConds(ii))=ii;
    end
    
    designMatrix     = conditions2design(designConditions);
    
    % Permute data and do denoising
    dataToDenoise = permute(data, [3,1,2]);
    [results, ~, ~, denoised_data] = ...
        denoisedata(designMatrix(~badEpochs,:),dataToDenoise(~badChannels, :, ~badEpochs), noisePoolFun, evokedfun, opt);
    
    % Add Nans back into denoised data
    dataDenoised = NaN(sz(3),sz(1),sz(2));
    dataDenoised(~badChannels, :, ~badEpochs) = denoised_data{1};
    
    % Plot figures
    if opt.verbose
        if ~exist(fullfile(dirPth.meg.saveFigPth,'postdenoise_spectra'), 'dir')
            mkdir(fullfile(dirPth.meg.saveFigPth,'postdenoise_spectra'));
        end
        
        figure;
        plot(results.origmodel.r2,results.finalmodel.r2,'k.'); hold on;
        plot([0 max([results.origmodel.r2,results.finalmodel.r2])], ...
            [0 max([results.origmodel.r2,results.finalmodel.r2])],'r-');
        ylabel('Final R2'); xlabel('Orig R2'); axis square;
        set(gca, 'TickDir', 'out', 'FontSize', 14);
        
        if opt.saveFig
            print(gcf, '-dpng', fullfile(dirPth.meg.saveFigPth, sprintf('%s_r2_denoising', subjID)))
        end
        
        freqIdx = mprfFreq2Index(size(denoised_data{1},2), opt.flickerFreq, opt.fs);
        t = (0:size(denoised_data{1},2)-1)./opt.fs;
        ft = (0:length(t)-1)/max(t);
        
        figure;
        for chan = 1:size(denoised_data{1},1)
            % Incoherent spectrum
            subplot(211); cla
            amps = abs(fft(denoised_data{1}(chan,:,:),[],2))/size(denoised_data{1},2)*2;
            
            plot(ft,nanmean(squeeze(amps),2)); hold on;
            plot([opt.flickerFreq opt.flickerFreq], [0 max(nanmean(squeeze(amps),2))]);
            title(sprintf('Incoherent spectrum: sensor %d', chan));
            xlim([1 100]); xlabel('Frequency (Hz)'); ylabel('Amplitudes (T)');
            set(gca, 'TickDir', 'out', 'FontSize', 14);
            
            drawnow;
            
            % Coherent spectrum
            subplot(212); cla
            meanDenoisedData = squeeze(nanmean(denoised_data{1}(chan,:,:),3));
            amps = abs(fft(meanDenoisedData,[],2))/size(meanDenoisedData,2)*2;
            
            plot(ft,amps); hold on;
            plot([opt.flickerFreq opt.flickerFreq], [min(amps) max(amps)]);
            title(sprintf('Coherent spectrum: sensor %d', chan));
            xlim([1 100]); xlabel('Frequency (Hz)'); ylabel('Amplitudes (T)');
            set(gca, 'TickDir', 'out', 'FontSize', 14);
            
            if opt.saveFig
                print(gcf, '-dpng', fullfile(dirPth.meg.saveFigPth,'postdenoise_spectra', ...
                    sprintf('%s_fft_spectrum_postDenoise_sensor%d', subjID,chan)))
            end
            drawnow;
        end
        
        % stim and blank SSVEF for incoherent spectrum
        allAmps                 = abs(fft(dataDenoised,[],2))/size(dataDenoised,2)*2;
        epochsStimToPlot        = triggers.stimConditions<10;
        stimDataToPlot.incoh    = squeeze(nanmean(allAmps(:,freqIdx,epochsStimToPlot),3));
        
        epochsBlankToPlot       = triggers.stimConditions==10;
        blankDataToPlot.incoh   = squeeze(nanmean(allAmps(:,freqIdx,epochsBlankToPlot),3));
        
        % stim and blank SSVEF for coherent spectrum
        meanStimTs  = nanmean(dataDenoised(:,:,epochsStimToPlot),3);
        meanBlankTs = nanmean(dataDenoised(:,:,epochsBlankToPlot),3);
        
        meanAmpsStim = abs(fft(meanStimTs, [], 2))/size(dataDenoised,2)*2;
        meanAmpsBlank = abs(fft(meanBlankTs, [], 2))/size(dataDenoised,2)*2;
        
        stimDataToPlot.coh = squeeze(meanAmpsStim(:,freqIdx));
        blankDataToPlot.coh = squeeze(meanAmpsBlank(:,freqIdx));

        fh1 = figure;
        megPlotMap(stimDataToPlot.coh,[],[],[],[],[],[],'interpmethod', 'nearest');
        title('Steady state visually evoked field (10 Hz); Coherent spectrum');
        
        fh2 = figure;
        megPlotMap(stimDataToPlot.coh-blankDataToPlot.coh,[],[],[],[],[],[],'interpmethod', 'nearest');
        title('Steady state visually evoked field - blanks (10 Hz) Coherent spectrum');
        
        fh3 = figure;
        megPlotMap(stimDataToPlot.incoh,[],[],[],[],[],[],'interpmethod', 'nearest');
        title('Steady state visually evoked field (10 Hz) Incoherent spectrum');
        
        fh4 = figure;
        megPlotMap(stimDataToPlot.incoh-blankDataToPlot.incoh,[],[],[],[],[],[],'interpmethod', 'nearest');
        title('Steady state visually evoked field - blanks (10 Hz) Incoherent spectrum');
        
        if opt.saveFig
            print(fh1, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_coh', subjID)))
            print(fh2, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_coh_diff', subjID)))
            print(fh3, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_incoh', subjID)))
            print(fh4, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_incoh_diff', subjID)))
        end
        
    end % opt.verbose
end % opt.doDenoise

%% 7. Reblock and save prepreocessed timeseries
% we reshape the data so that the matrix is will have 4 dimensions:
%   time x epochs x channels x runs,
% where epochs are again split by the number of runs

if opt.verbose; fprintf('(%s) Save data...\n', mfilename); end

if strcmp(subjID,'wlsubj030')
    assert(numRuns==10);
    assert(numOfEpochsPerRun==211);
else
    assert(numRuns==19);
    assert(numOfEpochsPerRun==140);
end

% Set up array for splitting up in blocks   (i.e. number of RUNS)
preprocDataBlocked =  NaN(sz(3), sz(1), sz(2)/numRuns, numRuns);

% plot all triggers
if opt.verbose
    figure;
    plot(triggers.stimConditions, 'LineWidth', 2);
    xlabel('Time (timepoints)'); ylabel('Trigger num'); hold all;
    if opt.saveFig
        print(gcf, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_triggers_stimCondition', subjID)))
    end
end

for n = 1:numRuns
    
    % Get first epoch
    startEpoch = startOfRun(n);
    
    % Get last epoch depending on the run number
    if n < numRuns
        lastEpoch  = startOfRun(n+1)-1;
        theseEpochs = startEpoch:lastEpoch;
    else
        lastEpoch  = size(preprocDataBlocked,3);
        theseEpochs = startEpoch:(startEpoch+lastEpoch-1);
    end
    
    % Get data and put in dataBlocked variable
    preprocDataBlocked(:,:,:,n) = dataDenoised(:,:, theseEpochs);
    
    % Mark those that are blocked
    if opt.verbose; plot(theseEpochs,triggers.stimConditions(theseEpochs),'r:', 'LineWidth', 4); end
end

% Clean up
clear data;

% Reshape
data.data = permute(preprocDataBlocked, [2, 3, 4, 1]); % channels x time x epochs x blocks --> time x epochs x blocks x channels

% Save if requested
if opt.doSaveData
    save(fullfile(saveDataPth, 'epoched_data_hp_preproc_denoised.mat'), 'data', '-v7.3');
end


cd(curDir);

end




