function [phRefAmp10Hz, bestRefPhase, maxVarExplVal, bestBetas, bestOffsets] = ...
    mprf_MEGPhaseReferenceData(megData, predMEGResponse, runGroup, opt, dirPth)
% Function to computing phase referenced amplitude from preprocessed MEG data
% and predicted MEG responses from cortical surface
%   phaseRefMEGResponse = mprf_MEGPhaseReferenceData(megData, predMEGResponse)
%
% There are two ways of computing the reference phase, depending on what
% FFT spectrum is used to extract the 10 Hz amplitudes and phases. If one
% uses the coherent spectrum (average data then computing fft) we use split-
% half-crossvalidation (by randomly dividing the runs into two groups).
% Data in halve is then used to compute one reference phase from the other.
% If we use the incoherent spectrum (compute fft then average data), we do
% a leave-one-out-crossvalidation, so use n-1 runs to compute the reference
% phase for the left out run.
%   We use different approaches depending on the spectral data, because the
% coherent spectrum needs to average a relatively large amount of data to
% remove noise and get a stable phase in either data halve, which is not
% the case for the incoherent spectrum.
%
% INPUTS:
%   megData         : preprocessed MEG data (time x epochs x run x sensors)
%   predMEGResponse : predicted MEG responses (epochs x sensors)
%   runGroup        : groups for split halves in case of coherent spectrum
%   opt             :  struct with boolean flag options

%
% OUTPUT:
%   phRefAmp10Hz        : Phase referenced MEG time series (epochs x runs x sensors)
%   bestRefPhase        : Ref phases that gives highest var explained (1 x runs x sensors)
%   maxVarExplVal       : Variance explained by best ref phases (1 x runs x sensors)
%   bestBetas           : corresponding beta (scale factor) of best predicting ref phase  (1 x runs x sensors)
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019


if ~exist('runGroup', 'var') || isempty(runGroup)
    runGroup = [];
end

% Define the number of references phases to try
phaseRange = linspace(0,2*pi,100); % range of values to search for the reference phase

% Check dimensions of MEG data
[nTimepoints, nEpochs, nRuns, nSensors] = size(megData);

% Make sure the first dimension is time (i.e. 1100 ms, the largest of all)
assert(all(nTimepoints>[nEpochs,nRuns,nSensors]))

% Get frequency index
freqIdx = mprfFreq2Index(nTimepoints, opt.meg.flickerFreq, opt.meg.fs);

% Check if we want the coherent (average before FFT) or incoherent spectrum
% (after after FFT)
if opt.meg.useCoherentSpectrum
    
    warning off
    fprintf('(%s): Checking best reference phase for coherent spectrum.\n', mfilename)
    fprintf('(%s): Using split half cross-validation..\n', mfilename)
    fprintf('(%s): Runs in split half A: %s\n', mfilename, string(num2str(runGroup{1})))
    fprintf('(%s): Runs in split half B: %s\n', mfilename, string(num2str(runGroup{2})))
    
    fprintf('Sensors:')
    for s = 1:nSensors
        
        fprintf('.%d',s)
        
        % Take the mean across runs first (time x epochs)
        meanTs{1} = nanmean(megData(:,:,runGroup{1},s),3);
        meanTs{2} = nanmean(megData(:,:,runGroup{2},s),3);
        
        % Transform to Fourier domain
        F{1} = fft(meanTs{1});
        F{2} = fft(meanTs{2});
        
        % Get phase and amplitudes
        ph{1}      = angle(F{1});
        amp{1}     = abs(F{1})/nTimepoints*2;
        
        ph{2}      = angle(F{2});
        amp{2}     = abs(F{2})/nTimepoints*2;
        
        % Select amplitude and phase at stimulus frequency (10 Hz)
        amp10Hz{1} = squeeze(amp{1}(freqIdx,:)); % one value per epoch
        ph10Hz{1}  = squeeze(ph{1}(freqIdx,:)); % one value per epoch
        
        amp10Hz{2} = squeeze(amp{2}(freqIdx,:)); % one value per epoch
        ph10Hz{2}  = squeeze(ph{2}(freqIdx,:)); % one value per epoch
        
        allAmp10Hz(:,:,s) = cat(2,amp10Hz{1}',amp10Hz{2}'); % one value per epoch x run x sensor
        allPh10Hz(:,:,s) = cat(2,ph10Hz{1}',ph10Hz{2}');
        
        currentnans{1} = isnan(amp10Hz{1});
        currentnans{2} = isnan(amp10Hz{2});
        
        for split = 1:length(runGroup)
            
            x0 = pi;
            A = amp10Hz{split}(~currentnans{split})';
            P = ph10Hz{split}(~currentnans{split})';
            pred = predMEGResponse(~currentnans{split},s);
            options = optimset('Display','off');
            
            if isempty(A)
                B = NaN; refph = NaN; offset = NaN;
            else
                refph = fminsearch(@(x) ...
                    phaseReferencedPrediction(x, pred, A, P, opt.addOffsetParam), x0,options);
                phRef10Hz = rescaleAmpsWithRefPhase(A, P, refph);
                [B, offset, ve] = regressPredictedResponse(phRef10Hz, pred, 'addOffsetParam', opt.addOffsetParam);
                if B < 0, B = -B; refph = refph + pi; end
            end
            
            bestBetas(1,split, s)     = B;
            bestRefPhase(1,split, s)  = mod(refph, 2*pi);
            bestOffsets(1,split, s)   = offset;
            maxVarExplVal(1,split, s) = ve;
            
        end % split halves
    end % sensors
    
    fprintf('...Done!\n')
    
else % if using incoherent spectrum (then start with FFT before averaging)
    
    % FFT of MEG data (freq, epochs, run, sensors)
    F = fft(megData);
    
    % Get phase and amplitudes
    ph      = angle(F);
    amp     = abs(F)/size(megData,1)*2;
    
    % Select amplitude and phase at stimulus frequency (10 Hz)
    allAmp10Hz = squeeze(amp(freqIdx,:,:,:)); % one value per epoch x run x sensor
    allPh10Hz  = squeeze(ph(freqIdx, :,:,:)); % one value per epoch x run x sensor
    
    varexpl           = NaN(length(phaseRange), nRuns, nSensors);
    refPhase          = NaN(length(phaseRange), nRuns, nSensors);
    
    warning off
    fprintf('(%s): Checking best reference phase for incoherent spectrum.\n', mfilename)
    fprintf('(%s): Using leave-one-out cross-validation.\n', mfilename)
    
    
    % Get test run (left-out run)
    for split = 1:nRuns
        fprintf('(%s): Runs in split half B: %s\n', mfilename, num2str(split))
        
        % Get leftout runs
        trainingRuns = setdiff(1:nRuns,split);
        fprintf('(%s): Training runs: %s\n', mfilename, string(num2str(trainingRuns)))
        
        fprintf('Sensors:')
        for s = 1:nSensors
            fprintf('.%d',s)
            
            currentnans = isnan(allAmp10Hz(:,trainingRuns,s))';
            
            % Select current ampl data from training runs
            A = mean(allAmp10Hz(~currentnans, trainingRuns, s),2); % (epochs x runs-1)
            
            % Get average phase across left out runs (epochs x 1)
            P = circularavg(allPh10Hz(~currentnans, trainingRuns, s), [], 2);
            
            x0 = pi;
            pred = predMEGResponse(~currentnans,s);
            options = optimset('Display','off');
            
            if isempty(A)
                B = NaN; refph = NaN; offset = NaN;
            else
                refph = fminsearch(@(x) ...
                    phaseReferencedPrediction(x, pred, A, P, opt.addOffsetParam), x0,options);
                phRef10Hz = rescaleAmpsWithRefPhase(A, P, refph);
                [B, offset, ve] = regressPredictedResponse(phRef10Hz, pred, 'addOffsetParam', opt.addOffsetParam);
                if B < 0, B = -B; refph = refph + pi; end
            end
            bestBetas(1,split, s)     = B;
            bestRefPhase(1,split, s)  = mod(refph, 2*pi);
            bestOffsets(1,split, s)   = offset;
            maxVarExplVal(1,split, s) = ve;
        end % sensors
        
    end % split halves
end % if opt.useCoherentSpectrum

%% Cross-validate phases,recompute phase-referenced amplitudes per sensor

% rescale the original amplitudes and phase from MEG data with best
% reference phase from other half or training data
if opt.meg.useCoherentSpectrum
    % make sure to swap the reference phases for half 1 and 2
    tmp = [bestRefPhase(1,2,:),bestRefPhase(1,1,:)];
    bestRefPhase = tmp;
    phRefAmp10Hz = rescaleAmpsWithRefPhase(allAmp10Hz, allPh10Hz, bestRefPhase);
else
    % not needed for the leave-one-out procedure, they are already
    % cross-validated with the leave-one-out-procedure
    phRefAmp10Hz = rescaleAmpsWithRefPhase(allAmp10Hz, allPh10Hz, bestRefPhase);
end
warning on

fprintf('(%s) Done!\n',mfilename)


if ~opt.vary.perturbOrigPRFs
    
    if opt.meg.doSplitHalfReliability % takes a long time, we default is set to false. See also 's_splithalfR_allSubjects.m'
        % Compute split half reliability of SSVEF amplitude halves
        splitHalfAmpCorrelation = mprf_ComputeSplitHalfReliability(dirPth, opt, megData, freqIdx);
    end
    
    if opt.verbose
        
        if opt.meg.doSplitHalfReliability
            % Plot split half amplitude reliability
            fH1 = figure; clf; megPlotMap(splitHalfAmpCorrelation,[0 max(splitHalfAmpCorrelation)],fH1, 'hot', ...
                'Mean split half reliability of SSVEF amplitudes', [],[], 'interpmethod', 'nearest');
            c = colorbar; c.Location='eastoutside';
            
            % Plot split half amplitude reliability
            fH2 = figure; clf; megPlotMap(splitHalfAmpCorrelation,[0 max(splitHalfAmpCorrelation)],fH2, 'hot', ...
                'Mean split half reliability of SSVEF amplitudes');
            c = colorbar; c.Location='eastoutside';
            
            if opt.saveFig
                figure(fH1)
                figurewrite(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'refphase', ...
                    'splitHalfAmplitudeReliability'),[],0,'.',1);
                figure(fH2)
                figurewrite(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'refphase', ...
                    'splitHalfAmplitudeReliability_interp'),[],0,'.',1);
            end
        end
        
        fH2 = figure(1); clf; set(gcf, 'Position',  [1000, 651, 1285, 687]);
        t = 1:size(phRefAmp10Hz,1);
        
        for s = 1:nSensors
            % Plot the split half amplitudes (absolute values, used in phase referencing)
            % do some plotting for debugging
            
            figure(fH2); clf;
            
            ax1 = subplot(2,1,1);
            plot(t, allAmp10Hz(:,1,s), 'r'); hold on; plot(t, allAmp10Hz(:,2,s), 'g');
            xlabel('time points'); ylabel('Magnetic flux (T)');
            title(sprintf('Sensor %d amplitudes',s)); %, splithalf reliability %1.2f', s, splitHalfAmpReliability(s)));
            legend({'Amplitudes of split half 1', 'Amplitudes of split half 2'}); box off;
            set(ax1, 'TickDir', 'out', 'FontSize', 10)
            
            % Plot the split half phase-referenced amplitudes
            ax2 = subplot(2, 1, 2);
            plot(t, phRefAmp10Hz(:,1,s), 'r'); hold on;
            plot(t, phRefAmp10Hz(:,2,s), 'g');
            plot(t, nanmean(phRefAmp10Hz(:,:,s),2), 'k:', 'lineWidth',3);
            title(sprintf('Best ref phases split halves: %1.2f %1.2f, resulting in %1.2f %1.2f var expl', bestRefPhase(:,:,s), maxVarExplVal(:,:,s)));
            plot(t, predMEGResponse(:,s).*bestBetas(:,1,s) + bestOffsets(:,1,s), 'b', 'lineWidth',3);
            plot(t, predMEGResponse(:,s).*bestBetas(:,2,s) + bestOffsets(:,2,s), 'c', 'lineWidth',3);
            xlabel('Time points'); ylabel('Magnetic flux (T)')
            legend({'Phase referenced split half 1', 'Phase referenced split half 2', ...
                'Phase ref mean', 'Predicted MEG resp (scaled with beta split half 1)' ...
                'Predicted MEG resp (scaled with beta split half 2)'}); box off;
            set(ax2, 'TickDir', 'out', 'FontSize', 10)
            
            if opt.saveFig
                % figurewrite(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'refphase', ...
                %     sprintf('sensor%d_amplitudes%s', s, opt.fNamePostFix)), [], [1 300], '.',1);
            end % savefig
        end % sensor
    end % verbose
end % non vary PRF analysis only

end


