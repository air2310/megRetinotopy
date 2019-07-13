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



% If perturb original pRFs, check dimensions with loaded pRF data
if strcmp(opt.perturbOrigPRFs, 'position')
    assert(size(predSurfaceResponse,3)==length(opt.varyPosition))
    nIter = length(opt.varyPosition);
    predSurfaceResponseAll = predSurfaceResponse; % Keep a copy of all responses
elseif strcmp(opt.perturbOrigPRFs, 'size')
    assert(size(predSurfaceResponse,3)==length(opt.varySize))
    nIter = length(opt.varySize);
    predSurfaceResponseAll = predSurfaceResponse; % Keep a copy of all responses
elseif strcmp(opt.perturbOrigPRFs, 'scramble')
    assert(size(predSurfaceResponse,3)==length(opt.nScrambles))
    nIter = length(opt.nScrambles);
    predSurfaceResponseAll = predSurfaceResponse; % Keep a copy of all responses
elseif ~opt.perturbOrigPRFs
    nIter = 1;
end


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
            print(fullfile(dirPth.model.saveFigPth, opt.subfolder, ...
                sprintf('predMEGResponseFromPRF%s_%d', opt.fNamePostFix, ii)), '-dpng')
        end
    end
end

% Remove last dimension out, if not used
predMEGResponse = squeeze(predMEGResponse);

if opt.doSaveData
    if ~exist(fullfile(dirPth.model.saveDataPth, opt.subfolder), 'dir')
        mkdir(fullfile(dirPth.model.saveDataPth, opt.subfolder)); end
    save(fullfile(dirPth.model.saveDataPth,'predMEGResponsesFromPRFs'),'predMEGResponse');
end

return
