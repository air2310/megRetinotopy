function phRefAmp10Hz = mprf_MEGPhaseReferenceData(megData, predMEGResponse, opt)
% Function to computing phase referenced amplitude from preprocessed MEG data
% and predicted MEG responses from cortical surface
%   phaseRefMEGResponse = mprf_MEGPhaseReferenceData(megData, predMEGResponse)
%
% INPUTS:
%   megData         : preprocessed MEG data (time x epochs x run x sensors)
%   predMEGResponse : predicted MEG responses (epochs x sensors)
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
fprintf('(%s) Checking best reference phase .', mfilename)
for rp = 1:length(phaseRange)
    fprintf('.')
    % Get reference phase
    thisRefPhase = phaseRange(rp);
    
    for run = 1:nRuns
        
        % Get leftout runs
        leftOutRuns = setdiff(1:nRuns,run);
        
        for sensor = 1:nSensors
            currentnans = isnan(amp10Hz(:,run,sensor));
            
            % Select current phase data from left out runs
            tmp = ph10Hz(:, leftOutRuns, sensor); % (epochs x runs-1)
            
            % Get average phase across left out runs (epochs x 1)
            leftOutAverage_ph10Hz = angle(nansum(exp(tmp*1i),2));
            
            % Compute difference between reference and actual run
            phaseDiff = leftOutAverage_ph10Hz - thisRefPhase;
            
            % Compute angle of phase difference
            angPhaseDiff = cos(phaseDiff);
            
            % Scale amplitude by angle of phase difference
            phRef10Hz = amp10Hz(~currentnans,run,sensor).*angPhaseDiff(~currentnans);
            
            % Get predicted MEG response and a set of ones for capturing
            % the average
            X = [ones(size(predMEGResponse(~currentnans,sensor))) predMEGResponse(~currentnans,sensor)];
            
            % Regress prediction from phase referenced 10 Hz MEG response
            [B,~,~,~,stats] = regress(phRef10Hz,X);
            
            varexpl(rp,run,sensor) = stats(1);
            
            if B(2) < 0
                % If regression results in a negative scale factor, then
                % add pi to the reference phase
                refPhase(rp,run,sensor) = thisRefPhase+pi;
            else
                refPhase(rp,run,sensor) = thisRefPhase;
            end
 
        end
    end
end
warning on
fprintf('\n(%s) done!\n',mfilename)

% Get phase that gives max CoD per run, per sensor
[~, maxVarExplIdx] = max(varexpl);

maxPhase     = refPhase(maxVarExplIdx);
phaseDiff    = ph10Hz - maxPhase;
maxAngle     = cos(phaseDiff);
phRefAmp10Hz = amp10Hz.*maxAngle;


return


