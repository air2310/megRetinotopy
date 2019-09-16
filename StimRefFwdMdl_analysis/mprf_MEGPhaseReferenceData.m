function [phRefAmp10Hz, bestRefPhase, maxVarExplVal, bestBetas] = mprf_MEGPhaseReferenceData(megData, predMEGResponse, runGroup, opt, dirPth)
% Function to computing phase referenced amplitude from preprocessed MEG data
% and predicted MEG responses from cortical surface
%   phaseRefMEGResponse = mprf_MEGPhaseReferenceData(megData, predMEGResponse)
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
%   betas
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
    fprintf('(%s): Checking best reference phase for coherent spectrum.', mfilename)
    
    for ll = 1:length(runGroup)
        fprintf('\nSplit half %d: ', ll)
        % Split all runs into halves
        leftInRuns = setdiff(1:nRuns,runGroup{ll}); % Left in --> data used to compute the best ref phase
        leftOutRuns = runGroup{ll}; % Left out --> data that will actually be rescaled with this computed ref phase
        
        fprintf('Sensors:')
        for s = 1:nSensors

            fprintf('.%d',s)

            % Take the mean across runs first (time x epochs)
            meanTs.in  = nanmean(megData(:,:,leftInRuns,s),3);
            meanTs.out = nanmean(megData(:,:,leftOutRuns,s),3);
            
            % Transform to Fourier domain
            F.in = fft(meanTs.in);
            F.out = fft(meanTs.out);
            
            % Get phase and amplitudes
            ph.in      = angle(F.in);
            amp.in     = abs(F.in)/nTimepoints*2;
            
            ph.out      = angle(F.out);
            amp.out     = abs(F.out)/nTimepoints*2;
            
            % Select amplitude and phase at stimulus frequency (10 Hz)
            amp10Hz.in = squeeze(amp.in(freqIdx,:)); % one value per epoch
            ph10Hz.in  = squeeze(ph.in(freqIdx, :)); % one value per epoch
            
            amp10Hz.out = squeeze(amp.out(freqIdx,:)); % one value per epoch
            ph10Hz.out  = squeeze(ph.out(freqIdx, :)); % one value per epoch
            
            allAmp10Hz(:,ll,s) = amp10Hz.out;
            allPh10Hz(:,ll,s) = ph10Hz.out;
            
            currentnans.in = isnan(amp10Hz.in);
            currentnans.out = isnan(amp10Hz.out);
            
            % Select current phase data from left out runs
            ph10Hz.in = ph10Hz.in(~currentnans.in); % (epochs x runs-10)
            
            for rp = 1:length(phaseRange)
                % Get reference phase
                thisRefPhase = phaseRange(rp);
                
                % Rescale amplitudes with diff between reference phase and
                % average phase of other runs
                phRef10Hz = rescaleAmpsWithRefPhase(amp10Hz.in(~currentnans.in), ph10Hz.in, thisRefPhase);
                
                % Regress prediction from phase referenced 10 Hz MEG response
                [B, ve] = regressPredictedResponse(phRef10Hz', predMEGResponse(~currentnans.in,s));
                betas(rp,ll,s) = B(2);
                varexpl(rp,ll,s) = ve;
                
                if B(2) < 0
                    % If regression results in a negative scale factor, then
                    % add pi to the reference phase
                    refPhase(rp,ll,s) = thisRefPhase+pi;
                else
                    refPhase(rp,ll,s) = thisRefPhase;
                end
                
            end % reference phase
        end % sensors
    end % left out runs
    
    
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
    fprintf('(%s): Checking best reference phase for incoherent spectrum .', mfilename)
    for rp = 1:length(phaseRange)
        fprintf('.')
        % Get reference phase
        thisRefPhase = phaseRange(rp);
        
        for leftOutRun = 1:nRuns
            
            % Get leftout runs
            otherRuns = setdiff(1:nRuns,leftOutRun);
            
            for sensor = 1:nSensors
                
                currentnans = isnan(allAmp10Hz(:,leftOutRun,sensor))';
                
                % Select current phase data from left out runs
                tmp = allPh10Hz(~currentnans, otherRuns, sensor); % (epochs x runs-1)
                
                % Get average phase across left out runs (epochs x 1)
                otherRunsAverage_ph10Hz = circularavg(tmp, [], 2);
                
                % Rescale amplitudes with diff between reference phase and
                % average phase of other runs
                phRef10Hz = rescaleAmpsWithRefPhase(allAmp10Hz(~currentnans,leftOutRun,sensor), otherRunsAverage_ph10Hz, thisRefPhase);
                
                % Regress prediction from phase referenced 10 Hz MEG response
                [B, ve] = regressPredictedResponse(phRef10Hz, predMEGResponse(~currentnans,sensor));
                
                betas(rp,ll,s) = B(2);
                varexpl(rp,leftOutRun,sensor) = ve;
                
                if B(2) < 0
                    % If regression results in a negative scale factor, then
                    % add pi to the reference phase
                    refPhase(rp,leftOutRun,sensor) = thisRefPhase+pi;
                else
                    refPhase(rp,leftOutRun,sensor) = thisRefPhase;
                end
                
            end % sensors
        end % left out runs
    end % reference phase
end  % if opt.useCoherentSpectrum

%% Get phase that gives max variance explained per run, per sensor

% Allocate space
bestBetas = NaN(1,nSensors,size(varexpl,2));
bestRefPhase = NaN(1,nSensors,size(varexpl,2));

% Loop over run (group) to make sure the indexing will be correct
for rr = 1:size(varexpl,2)
    
    varexplRun = squeeze(varexpl(:,rr,:));
    betaRun    = squeeze(betas(:,rr,:));
    refPhaseRun = squeeze(refPhase(:,rr,:));
    
    % Get max variance explained
    [maxVEvalue, maxVEidx]  = max(varexplRun, [], 'omitnan');
    
    % Get nan sensors
    nanIdx = isnan(maxVEvalue);
    
    maxVarExplIdx(rr,:,:) = maxVEidx;
    maxVarExplVal(rr,:,:) = maxVEvalue;

    bestBetas(1,~nanIdx,rr)     = betaRun(maxVEidx(~nanIdx))';
    bestRefPhase(1,~nanIdx,rr)  = refPhaseRun(maxVEidx(~nanIdx))';
    
end

maxVarExplIdx = permute(maxVarExplIdx, [2,1,3]);
maxVarExplVal = permute(maxVarExplVal, [2,1,3]);
bestBetas     = permute(bestBetas, [1,3,2]);
bestRefPhase  = permute(bestRefPhase, [1,3,2]);

% rescale the original amplitudes and phase from MEG data
phRefAmp10Hz = rescaleAmpsWithRefPhase(allAmp10Hz, allPh10Hz, bestRefPhase);


warning on
fprintf('\n(%s) done!\n',mfilename)

% do some plotting for debugging
if ~opt.vary.perturbOrigPRFs
    fH1 = figure(1); set(gcf, 'Position',  [1000, 651, 1285, 687]);
    
    for s = 1:nSensors
        figure(fH1); clf;
        
        ax1 = subplot(2,1,1);
        plot(1:140, allAmp10Hz(:,1,s), 'r'); hold on; plot(1:140, allAmp10Hz(:,2,s), 'g');
        xlabel('time points'); ylabel('Magnetic flux (T)'); title(sprintf('Sensor %d amplitudes', s));
        legend({'Amplitudes of split half 1', 'Amplitudes of split half 2'}); box off;
        set(ax1, 'TickDir', 'out', 'FontSize', 10)
        
        ax2 = subplot(2, 1, 2);
        plot(1:140, phRefAmp10Hz(:,1,s), 'r'); hold on; 
        plot(1:140, phRefAmp10Hz(:,2,s), 'g');
        plot(1:140, nanmean(phRefAmp10Hz(:,:,s),2), 'k:', 'lineWidth',3);
        title(sprintf('Best ref phases split halves: %1.2f %1.2f, resulting in %1.2f %1.2f var expl', bestRefPhase(:,:,s), maxVarExplVal(:,:,s)));
        plot(1:140, predMEGResponse(:,s).*bestBetas(:,2,s), 'b');
        xlabel('time points'); ylabel('Magnetic flux (T)')
        legend({'Phase referenced split half 1', 'Phase referenced split half 2', ...
            'Phase ref mean', 'Predicted MEG resp (scaled with beta)'}); box off;
        set(ax2, 'TickDir', 'out', 'FontSize', 10)
        
        print(fH1,fullfile(dirPth.model.saveFigPth, opt.subfolder, 'refphase', ...
            sprintf('sensor%d_amplitudes%s', s, opt.fNamePostFix)), '-dpng')
        
    end
end

return


