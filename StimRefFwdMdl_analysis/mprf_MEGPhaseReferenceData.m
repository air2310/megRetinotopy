function phaseRefMEGResponse = mprf_MEGPhaseReferenceData(megData, predMEGResponse, opt)
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
phaseRange = linspace(0,pi,20); % range of values to search for the reference phase

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

phRef10Hz = NaN(length(phaseRange), nEpochs, nRuns, nSensors);
B         = NaN(length(phaseRange), nRuns, nSensors);
predPhRefResponse = NaN(length(phaseRange), nEpochs, nRuns, nSensors);
CoD       = NaN(length(phaseRange), nRuns, nSensors);

for rp = 1:length(phaseRange)
    
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
            phRef10Hz(rp, ~currentnans,run,sensor) = amp10Hz(~currentnans,run,sensor).*angPhaseDiff(~currentnans);
            
            % Get predicted MEG response and a set of ones for capturing
            % the average
            X = [ones(size(predMEGResponse(~currentnans,sensor))) predMEGResponse(~currentnans,sensor)];
            
            % Regress prediction from phase referenced 10 Hz MEG response
            B_tmp = X \ phRef10Hz(rp, ~currentnans,run,sensor)';
            
            % Pick reference phase that results in positive amplitudes, and thus
            % positive beta values for stimulus peaks in predicted MEG response
            if B_tmp(2)<0
                fprintf('(%s): RefPhase %d, Sensor %d, Run %d, current beta is negative: %d \n', mfilename, thisRefPhase, sensor, run, B_tmp(2)) % skip
            else
                B(rp,run,sensor) = B_tmp(2);
            
                % Get the predicted times series with this reference phase beta
                predPhRefResponse(rp,~currentnans,run,sensor) =  X * B_tmp;

                % Compute coefficient of determination:
                CoD(rp,run,sensor) = 1 - (var(phRef10Hz(rp, ~currentnans,run,sensor) - predPhRefResponse(rp,~currentnans,run,sensor)) ...
                                ./ var(phRef10Hz(rp, ~currentnans,run,sensor)));
            end
        end
    end
end

% Get phase that gives max CoD per run, per sensor
[maxPhase, maxPhaseIdx] = max(CoD);

phaseDiff    = ph10Hz - maxPhase;
maxAngle     = cos(phaseDiff);
phRefAmp10Hz = amp10Hz.*maxAngle;

% Take mean across 19 runs
meanPhRefAmp10Hz = squeeze(nanmean(phRefAmp10Hz,2));


%% Get predicted response that explains most variance from mean response

% preallocate space
meanPredResponse = NaN(nEpochs,nSensors);
meanCoD = NaN(1,nSensors);
for s = 1:nSensors
    
    % Identify and remove nans
    meanNanMask = isnan(meanPhRefAmp10Hz(:,s));
    meanPrediction = predMEGResponse(~meanNanMask,s);
    meanData = meanPhRefAmp10Hz(~meanNanMask,s);
    
    % Create predictions
    meanX = [ones(size(meanPrediction)) meanPrediction];
    
    % Regress out predictions
    meanB = meanX \ meanData;
    
    % Compute scaled predictions with betas
    meanPredResponse(~meanNanMask,s) =  meanX * meanB;

    % Compute coefficient of determination:
    meanCoD(s) = 1 - (var(meanData - meanPredResponse(~meanNanMask,s)) ...
                                ./ var(meanData));

end

return


