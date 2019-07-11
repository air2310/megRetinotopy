function phRefAmp10Hz = mprf_MEGPhaseReferenceData(megData, predMEGResponse, dirPth, opt)
% Function to computing phase referenced amplitude from preprocessed MEG data
% and predicted MEG responses from cortical surface
%   phaseRefMEGResponse = mprf_MEGPhaseReferenceData(megData, predMEGResponse)
%
% INPUTS:
%   megData         : preprocessed MEG data (time x epochs x run x sensors)
%   predMEGResponse : predicted MEG responses (epochs x sensors)
%   dirPth          : paths to files for given subject
%   opt             :  struct with boolean flag options

%
% OUTPUT:
%   phaseRefMEGResponse : Phase referenced MEG time series (sensors x epochs)


% Define the number of references phases to try
phaseRange = linspace(0,2*pi,100); % range of values to search for the reference phase

% Check dimensions of MEG data
[nTimepoints, nEpochs, nRuns, nSensors] = size(megData);


% Make sure the first dimension is time (i.e. 1100 ms, the largest of all)
assert(all(nTimepoints>[nEpochs,nRuns,nSensors]))

% FFT of MEG data (freq, epochs, run, sensors)
F = fft(megData);

% Get phase and amplitude information,
freqIdx = mprfFreq2Index(nTimepoints, opt.flickerFreq, opt.fs);

ph      = angle(F);
amp     = abs(F)/size(megData,1)*2;

% Select amplitude and phase at stimulus frequency (10 Hz)
amp10Hz = squeeze(amp(freqIdx,:,:,:)); % one value per epoch x run x sensor
ph10Hz  = squeeze(ph(freqIdx, :,:,:)); % one value per epoch x run x sensor

varexpl           = NaN(length(phaseRange), nRuns, nSensors);
refPhase          = NaN(length(phaseRange), nRuns, nSensors);

warning off
fprintf('(%s): Checking best reference phase .', mfilename)
for rp = 1:length(phaseRange)
    fprintf('.')
    % Get reference phase
    thisRefPhase = phaseRange(rp);
    
    for leftOutRun = 1:nRuns
        
        % Get leftout runs
        otherRuns = setdiff(1:nRuns,leftOutRun);
        
        for sensor = 1:nSensors
            
            currentnans = isnan(amp10Hz(:,leftOutRun,sensor))';
            
            % Select current phase data from left out runs
            tmp = ph10Hz(~currentnans, otherRuns, sensor); % (epochs x runs-1)
            
            % Get average phase across left out runs (epochs x 1)
            otherRunsAverage_ph10Hz = circularavg(tmp, [], 2);
            
            % Rescale amplitudes with diff between reference phase and
            % average phase of other runs
            phRef10Hz = rescaleAmpsWithRefPhase(amp10Hz(~currentnans,leftOutRun,sensor), otherRunsAverage_ph10Hz, thisRefPhase);
            
             % Regress prediction from phase referenced 10 Hz MEG response
            [B, ve] = regressPredictedResponse(phRef10Hz, predMEGResponse(~currentnans,sensor));
            
            varexpl(rp,leftOutRun,sensor) = ve;
            
            if B(2) < 0
                % If regression results in a negative scale factor, then
                % add pi to the reference phase
                refPhase(rp,leftOutRun,sensor) = thisRefPhase+pi;
            else
                refPhase(rp,leftOutRun,sensor) = thisRefPhase;
            end
            
        end
    end
end
warning on
fprintf('\n(%s) done!\n',mfilename)

% Get phase that gives max CoD per run, per sensor
[maxVarExplVal, maxVarExplIdx] = max(varexpl);

bestRefPhase = refPhase(maxVarExplIdx);
phRefAmp10Hz = rescaleAmpsWithRefPhase(amp10Hz, ph10Hz, bestRefPhase);


if opt.verbose
    
    % Visualize the mean reference phase across sensors
    fH1 = figure; clf;
    megPlotMap(circularavg(squeeze(bestRefPhase),[],1), [0 2*pi],[], 'hsv','Mean Ref phase across 19 runs of MEG data', [],[],'interpmethod', 'nearest')
    
    % Visualize the variance explained for every run across sensors
    fH2 = figure; set(fH2, 'Position', [17,578,1543,760]);
    for r = 1:nRuns
        subplot(3,7,r);
        dataToPlot = squeeze(maxVarExplVal(:,r,:));
        if all(isnan(dataToPlot))
            dataToPlot = zeros(size(dataToPlot));
        end
        megPlotMap(dataToPlot, [0 0.6],[], 'parula',sprintf('Run %d', r));
        drawnow;
    end
    
    if opt.saveFig
        print(fH1,fullfile(dirPth.model.saveFigPth, ...
            sprintf('mesh_xvalRefPhase_benson%d_highres%d_smoothed%d', ...
            opt.useBensonMaps, opt.fullSizeGainMtx, opt.useSmoothedData)), '-dpng')
        
        print(fH2,fullfile(dirPth.model.saveFigPth, ...
            sprintf('mesh_varExplbyXvalRefPhase_benson%d_highres%d_smoothed%d', ...
            opt.useBensonMaps, opt.fullSizeGainMtx, opt.useSmoothedData)), '-dpng')
    end
    
    % Visualize the reference phase for every run across sensors
    fH3 = figure;
    for s = 1:nSensors
        clf;
        mprf_polarplot(ones(size(bestRefPhase,2),1),bestRefPhase(1,:,s));
        title(sprintf('Best xval reference phases for 19 runs, sensor %d - mean r^2: %1.2f',s,nanmean(maxVarExplVal(1,:,s),2)));
        drawnow;
        if opt.saveFig
            if ~exist(fullfile(dirPth.model.saveFigPth, 'refphase'),'dir')
                mkdir(fullfile(dirPth.model.saveFigPth, 'refphase')); end
            print(fH3,fullfile(dirPth.model.saveFigPth, 'refphase', ...
                sprintf('sensor%d_xvalRefPhase_benson%d_highres%d_smoothed%d', ...
                s, opt.useBensonMaps, opt.fullSizeGainMtx, opt.useSmoothedData)), '-dpng')
        end
    end
end % opt.verbose

return


