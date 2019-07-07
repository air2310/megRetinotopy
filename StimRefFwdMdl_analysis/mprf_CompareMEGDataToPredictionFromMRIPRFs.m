function [meanPredResponse,meanVarExpl] = mprf_CompareMEGDataToPredictionFromMRIPRFs(phRefAmp10Hz, predMEGResponse, dirPth, opt)
% Function to compare phase referenced steady-state data MEG data
% to predicted MEG responses from MRI prfs
%
%    [meanPredResponse,meanVarExpl] = mprf_CompareMEGDataToPredictionFromMRIPRFs(phRefAmp10Hz, predMEGResponse, opt)
%
% INPUTS:
%   phRefAmp10Hz    : phase-referenced steady-state MEG sensor data (epochs x run x sensors)
%   predMEGResponse : predicted MEG responses (epochs x sensors)
%   opt             : struct with options
%
% OUTPUT:
%   meanPredResponse : mean predicted response (sensors x epochs)
%   meanVarExpl      : variance explained of mean data by modelfit (1 x sensors)


% Take mean across 19 runs
meanPhRefAmp10Hz = squeeze(nanmean(phRefAmp10Hz,2));

[nEpochs, ~, nSensors] = size(phRefAmp10Hz);

% Preallocate space
meanPredResponse = NaN(nEpochs,nSensors);
meanVarExpl      = NaN(1,nSensors);

% Get predicted response that explains most variance from mean response
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
    meanVarExpl(s) = 1 - (var(meanData - meanPredResponse(~meanNanMask,s)) ...
        ./ var(meanData));
    
end

% Plot some figures
if opt.verbose
    
    % remove nan sensors and sort by var expl.
    [val, idx] = sort(meanVarExpl, 'descend');
    tmp = idx(~isnan(val)); top10=tmp(1:10);
    ve = val(~isnan(val));
    
    % Plot var expl mesh
    figure; megPlotMap(meanVarExpl,[0 max(meanVarExpl)],[], 'parula','Var expl of mean phase-ref MEG data by modelfit');
    if opt.saveFig 
        if ~exist(dirPth.model.saveFigPth); mkdir(dirPth.model.saveFigPth); end
        print(fullfile(dirPth.model.saveFigPth, sprintf('varexpl_mesh_highres%d_smoothed%d_bensonmaps%d',opt.fullSizeGainMtx, opt.useSmoothedData, opt.useBensonMaps)), '-dpng');
    end
    
    % Plot Mean phase-referenced steady-state response and predicted response to
    % stimulus for top 10 sensors
    t = (0:nEpochs-1) .* diff(opt.epochStartEnd);
    
    figure; clf; set(gcf, 'Position', [652, 38,1206,1300], 'Color', 'w', 'Name', ...
        'Mean phase-ref MEG data and predicted response from pRF');
    
    for ii = 1:length(top10)
        subplot(5,2,ii);
        plot(t, meanPhRefAmp10Hz(:,top10(ii)), 'ko-', 'LineWidth',2);
        hold on; plot(t, meanPredResponse(:,top10(ii)), 'r', 'LineWidth',4);
        title(sprintf('Sensor %d, var expl: %1.2f',top10(ii), ve(ii)))
        xlabel('Time (s)'); ylabel('MEG response (Tesla)');
        set(gca, 'FontSize', 14, 'TickDir','out'); box off
    end
    
    if opt.saveFig
        print(fullfile(dirPth.model.saveFigPth, sprintf('varexpl_timeseries_highres%d_smoothed%d_bensonmaps%d',opt.fullSizeGainMtx, opt.useSmoothedData, opt.useBensonMaps)), '-dpng');
    end
    
end


return