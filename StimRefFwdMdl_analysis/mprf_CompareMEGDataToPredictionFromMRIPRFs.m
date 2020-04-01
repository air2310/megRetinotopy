function [predMEGResponseScaled, meanVarExpl, meanPhRefAmp10Hz] = ...
    mprf_CompareMEGDataToPredictionFromMRIPRFs(phRefAmp10Hz, predMEGResponse, bestBetas, bestOffsets, runGroup, opt)

% Function to compare phase referenced steady-state data MEG data
% to predicted MEG responses from MRI prfs
%
%    [meanPredResponse,meanVarExpl] = mprf_CompareMEGDataToPredictionFromMRIPRFs(phRefAmp10Hz, predMEGResponse)
%
% INPUTS:
%   phRefAmp10Hz    : phase-referenced steady-state MEG sensor data (epochs x run x sensors)
%   predMEGResponse : predicted MEG responses (epochs x sensors)
%
% OUTPUTS:
%   predMEGResponseScaled : mean predicted response scaled by beta (+offset) (epochs x sensors)
%   meanVarExpl           : variance explained of mean data by modelfit (1 x sensors)
%   meanPhRefAmp10Hz      : mean phase-referenced steady-state MEG sensors response
%                           (epochs x sensors)
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019


% Take mean across 19 runs
meanPhRefAmp10Hz = squeeze(nanmean(phRefAmp10Hz,2));

% Check dimensions
[~, nSensors] = size(meanPhRefAmp10Hz);

% Preallocate space
meanVarExpl           = NaN(1,nSensors);
predMEGResponseScaled = NaN(size(predMEGResponse));

% Get predicted response that explains most variance from mean response
if opt.recomputeFinalPredictionBetas
    for s = 1:nSensors
        [beta, offset, ~] = regressPredictedResponse(meanPhRefAmp10Hz(:,s), predMEGResponse(:,s), 'regressionType', opt.regressionType);
        
        predMEGResponseScaled(:,s) = predMEGResponse(:,s) .*beta + offset;
        
        % Compute coefficient of determination:
        meanVarExpl(s) = 1 - (  sum( (meanPhRefAmp10Hz(:,s) - predMEGResponseScaled(:,s)).^2, 'omitnan') ...
            ./ sum(meanPhRefAmp10Hz(:,s).^2, 'omitnan'));
    end
    
    
else
    % Select prf parameters of current iteration, scale with appropiate
    % betas and then average across all runs
    nRuns = length(runGroup{1})+length(runGroup{2});
    weights = [length(runGroup{1}) length(runGroup{2})]./nRuns;
    for s = 1:nSensors

        predMEGResponseScaled1 = predMEGResponse(:,s) .*bestBetas(1,s) + bestOffsets(1,s);
        predMEGResponseScaled2 = predMEGResponse(:,s) .*bestBetas(2,s) + bestOffsets(2,s);
    
        predMEGResponseScaled(:,s) = nansum(cat(3,weights(1).*predMEGResponseScaled1,weights(2).*predMEGResponseScaled2),3);
    
        % Compute coefficient of determination:
        meanVarExpl(s) = 1 - (  sum( (meanPhRefAmp10Hz(:,s) - predMEGResponseScaled(:,s)).^2, 'omitnan') ...
            ./ sum(meanPhRefAmp10Hz(:,s).^2, 'omitnan'));
    end
    
end



return