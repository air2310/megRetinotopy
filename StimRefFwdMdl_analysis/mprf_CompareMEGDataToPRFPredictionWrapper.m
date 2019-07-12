function [meanPredResponse,meanVarExpl] = mprf_CompareMEGDataToPRFPredictionWrapper(phRefAmp10Hz, predMEGResponse, dirPth, opt)
% Function to compare phase referenced steady-state data MEG data to
% predicted MEG responses from MRI prfs
%
%    [meanPredResponse,meanVarExpl] =
%    mprf_CompareMEGDataToPRFPredictionWrapper(phRefAmp10Hz, ...
%    predMEGResponse, opt)
%
% INPUTS:
%   phRefAmp10Hz     : phase-referenced steady-state MEG sensor data
%                       (epochs x run x sensors)
%   predMEGResponse  : predicted MEG responses
%                       (epochs x sensors)
%   dirPth           : paths to files for given subject
%   opt              : struct with options
%
% OUTPUT:
%   meanPredResponse : mean predicted response
%                       (sensors x epochs x optional variations)
%   meanVarExpl      : variance explained of mean data by modelfit 
%                       (1 x sensors [x optional variations])
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019



% If perturb original pRFs, check dimensions with loaded pRF data
if strcmp(opt.perturbOrigPRFs, 'position')
    assert(size(phRefAmp10Hz,4)==length(opt.varyPosition))
    assert(size(predMEGResponse,3)==length(opt.varyPosition))
    nIter = length(opt.varyPosition);
elseif strcmp(opt.perturbOrigPRFs, 'size')
    assert(size(phRefAmp10Hz,4)==length(opt.varySize))
    assert(size(predMEGResponse,3)==length(opt.varySize))
    nIter = length(opt.varySize);
elseif strcmp(opt.perturbOrigPRFs, 'scramble')
    assert(size(phRefAmp10Hz,4)==length(opt.nScrambles))
    assert(size(predMEGResponse,3)==length(opt.nScrambles))
    nIter = length(opt.nScrambles);
elseif ~opt.perturbOrigPRFs
    nIter = 1;
end

% Keep a copy of all responses
if opt.perturbOrigPRFs
    predMEGResponseAll = predMEGResponse;
    phRefAmp10HzAll    = phRefAmp10Hz;
end

% Allocate space
[~, nEpochs, ~, nSensors] = size(megData);

meanPredResponse = NaN(nSensors, nEpochs, nIter);
meanVarExpl      = NaN(nSensors, nIter);


% loop over dimensions, if necessary
for ii = 1:nIter
    
    % Select new prf parameters, if they vary in size or position
    predMEGResponse = predMEGResponseAll(:,:,ii);
    phRefAmp10Hz    = phRefAmp10HzAll(:,:,:,ii);
    
    [meanPredResponse(:,:,ii),meanVarExpl(:,ii)] = ...
        mprf_CompareMEGDataToPredictionFromMRIPRFs(phRefAmp10Hz, predMEGResponse);
    
    
    
    
    %% Plot summary figures
    if opt.verbose
        
        % remove nan sensors and sort by var expl.
        [val, idx] = sort(meanVarExpl(:,ii), 'descend');
        tmp = idx(~isnan(val)); top10=tmp(1:10);
        ve = val(~isnan(val));
        
        ttlPostFix = strsplit(sprintf('%s',opt.fNamePostFix), '_');
        ttl = sprintf('Var expl of mean phase-ref MEG data by modelfit %d %s %s %s', ii, ttlPostFix{2},ttlPostFix{3},ttlPostFix{4});
        
        % Plot var expl mesh
        figure; megPlotMap(meanVarExpl,[0 max(meanVarExpl)],[], 'parula', ttl, [],[], 'interpmethod', 'nearest');
        if opt.saveFig
            print(fullfile(dirPth.model.saveFigPth, opt.subfolder, sprintf('varexpl_mesh%s_%d',opt.fNamePostFix, ii)), '-dpng');
        end
        
        % Plot Mean phase-referenced steady-state response and predicted response to
        % stimulus for top 10 sensors
        t = (0:nEpochs-1) .* diff(opt.epochStartEnd);
        
        figure; clf; set(gcf, 'Position', [652, 38,1206,1300], 'Color', 'w', 'Name', ...
            sprintf('Mean phase-ref MEG data and predicted response from pRF %d %s %s %s', ii, ttlPostFix{2},ttlPostFix{3},ttlPostFix{4}));
        
        for tt = 1:length(top10)
            subplot(5,2,tt);
            plot(t, meanPhRefAmp10Hz(:,top10(tt)), 'ko-', 'LineWidth',2);
            hold on; plot(t, meanPredResponse(:,top10(tt)), 'r', 'LineWidth',4);
            title(sprintf('Sensor %d, var expl: %1.2f',top10(tt), ve(tt)))
            xlabel('Time (s)'); ylabel('MEG response (Tesla)');
            set(gca, 'FontSize', 14, 'TickDir','out'); box off
            ylim([-6,6].*10^-14); xlim([0, length(t)])
            legend({'Data', 'Prediction'}, 'Location', 'SouthWest'); legend boxoff
        end
        
        if opt.saveFig
            print(fullfile(dirPth.model.saveFigPth, opt.subfolder, sprintf('varexpl_timeseries_TOP10%s_%d',opt.fNamePostFix, ii)), '-dpng');
        end
        
        
        % Plot all timeseries separately
        figure; set(gcf, 'Position', [652,938,884,400], 'Color', 'w', 'Name', ...
            sprintf('Mean phase-ref MEG data and predicted response from pRF %d %s %s %s', ii, ttlPostFix{2},ttlPostFix{3},ttlPostFix{4}));
        
        for s = 1:nSensors
            clf;
            plot(t, meanPhRefAmp10Hz(:,s), 'ko-', 'LineWidth',2);
            hold on; plot(t, meanPredResponse(:,s), 'r', 'LineWidth',4);
            title(sprintf('Sensor %d, var expl: %1.2f',s, meanVarExpl(s)))
            xlabel('Time (s)'); ylabel('MEG response (Tesla)');
            set(gca, 'FontSize', 14, 'TickDir','out'); box off
            ylim([-6,6].*10^-14); xlim([0, length(t)])
            legend({'Data', 'Prediction'}, 'Location', 'SouthWest'); legend boxoff
            if opt.saveFig
                if ~exist(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'timeseries'), 'dir') 
                    mkdir(fullfile(dirPth.model.saveFigPth, subfolder, 'timeseries')); end
                print(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'timeseries', sprintf('varexpl_timeseries_sensor%d%s_iter%d',s, opt.fNamePostFix,ii)), '-dpng');
            end
        end
    end
end