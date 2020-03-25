function [meanVarExpl, meanPhRefAmp10Hz] = mprf_CompareMEGDataToPredictionFromMRIPRFs(phRefAmp10Hz, predMEGResponse)
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
%   meanPredResponse : mean predicted response (epochs x sensors)
%   meanVarExpl      : variance explained of mean data by modelfit (1 x sensors)
%   meanPhRefAmp10Hz : mean phase-referenced steady-state MEG sensors response
%                        (epochs x sensors)
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019


% Take mean across 19 runs
meanPhRefAmp10Hz = squeeze(nanmean(phRefAmp10Hz,2));   

% Check dimensions
[~, nSensors] = size(meanPhRefAmp10Hz);

% Preallocate space
meanVarExpl      = NaN(1,nSensors);

% Get predicted response that explains most variance from mean response
for s = 1:nSensors
    
    % Compute coefficient of determination:
    meanVarExpl(s) = 1 - (var(meanPhRefAmp10Hz(:,s) - predMEGResponse(:,s), 'omitnan') ...
        ./ var(meanPhRefAmp10Hz(:,s), 'omitnan')); 
end



return