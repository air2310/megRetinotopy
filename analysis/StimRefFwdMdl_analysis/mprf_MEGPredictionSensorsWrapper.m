function predMEGResponse = mprf_MEGPredictionSensorsWrapper(predSurfaceResponse, gain, dirPth, opt)
% Wrapper function for mprf_MEGPredictionSensors. This wrapper was
% created to keep the function to make predictions from prf data simple,
% while still allowing multiple analysis options.
%
% INPUTS:
%   predSurfaceResponse : predicted responses from surface
%                           (epochs x vertices)
%   gain                : gain matrix, defines weighted sum of vertices
%                           contributing to each individual MEG sensor
%                           (sensors x vertices)
%   dirPth              : paths to files for given subject
%   opt                 : struct with boolean flag options
%
% OUTPUT:
%   predMEGResponse   : predicted response time series for every MEG sensor
%                          (epochs x sensors x optional variations)
%
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019

% Check dimensions with loaded pRF data, and set the number of iterations
iter  = checkNumberOfIterations(predSurfaceResponse, opt, 'prfSensor');
nIter = length(iter);

% Keep a copy of all responses
predSurfaceResponseAll = predSurfaceResponse;

% Allocate space
predMEGResponse = NaN(size(predSurfaceResponse,1), size(gain,1) ,nIter);

% loop over dimensions, if necessary
for ii = 1:nIter
    
    % Select new prf parameters, if they vary in size or position
    predSurfaceResponse = predSurfaceResponseAll(:,:,ii);
    
    % Get predicted MEG response for this iteration
    predMEGResponse(:,:,ii) = mprf_MEGPredictionSensors(predSurfaceResponse, gain);
end

%% Debug figures
if opt.verbose
    if opt.saveFig
        if ~exist(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'predMEGResponse'), 'dir')
            mkdir(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'predMEGResponse'));
        end
    end
    
    figure, set(gcf, 'Position', [652   784   908   554], 'Color', 'w');
    for ii = 1:nIter
        clf;
        plot(predMEGResponse(:,:,ii));
        title(sprintf('Predicted MEG response from MRI prfs %d', ii))
        xlabel('Time points (epochs)');
        ylabel('MEG response (T)');
        ylim([-1,1]*max(abs(predMEGResponse(:))))
        set(gca, 'FontSize', 14, 'TickDir', 'out'); box off;
        
        if opt.saveFig
            print(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'predMEGResponse', ...
                sprintf('predMEGResponseFromPRF%s_%d', opt.fNamePostFix, ii)), '-dpng')
        end
    end
    
    figure, set(gcf, 'Position', [652   784   908   554], 'Color', 'w');
    if nIter > 1
        cols = round(nIter/3);
        rows = round(nIter/cols)+1;
    else
        cols = 1; rows = 1;
    end
    for ii = 1:nIter
        subplot(cols,rows,ii)
        megPlotMap(squeeze(nanvar(predMEGResponse(:,:,ii))),[],[],[],[],[],[],'interpmethod', 'nearest');
    end
    if opt.saveFig
        print(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'predMEGResponse', ...
            sprintf('predMEGResponseFromPRF%s_varianceAcrossEpochs', opt.fNamePostFix)), '-dpng')
    end
    
end

% Remove last dimension out, if not used
predMEGResponse = squeeze(predMEGResponse);

if opt.saveData
    if ~exist(fullfile(dirPth.model.saveDataPth, opt.subfolder, 'pred_resp'), 'dir')
        mkdir(fullfile(dirPth.model.saveDataPth, opt.subfolder, 'pred_resp')); end
    save(fullfile(dirPth.model.saveDataPth,opt.subfolder, 'pred_resp','predMEGResponsesFromPRFs'),'predMEGResponse', '-v7.3');
end

return
