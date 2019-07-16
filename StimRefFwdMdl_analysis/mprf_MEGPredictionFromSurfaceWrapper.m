function predSurfResponse = mprf_MEGPredictionFromSurfaceWrapper(prfSurfPath, stim, dirPth, opt)
% Wrapper function for mprf_MEGPredictionFromSurface. This wrapper was
% created to keep the function to make predictions from prf data simple,
% while still allowing multiple analysis options.
%
% INPUTS:
%   prfSurfPath     : directory where prf data live (string)
%   stim            : meg stimulus struct (should contain fieldname 'im', with windowed stimulus images for every epoch)
%   dirPth          : paths to files for given subject (struct)
%   opt             : struct with boolean flags
%
% OUTPUTS:
%   predSurfResponse : predicted response time series (vertex x timepoints [x optional variations])
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019



% Create folders if figures are saved
if opt.saveFig
    if ~exist(fullfile(dirPth.model.saveFigPth, opt.subfolder),'dir')
        mkdir(fullfile(dirPth.model.saveFigPth, opt.subfolder));
    end
end


% Define pRF parameters to read from surface
if strcmp(opt.perturbOrigPRFs, 'position')
    prfParams = {'varexplained', 'mask', 'recomp_beta', 'x_vary.mgz', 'y_vary.mgz', 'sigma_smoothed'};
elseif strcmp(opt.perturbOrigPRFs, 'size')
    prfParams = {'varexplained', 'mask', 'recomp_beta', 'x_smoothed', 'y_smoothed', 'sigma_vary.mgz'};
elseif strcmp(opt.perturbOrigPRFs, 'scramble')
    prfParams = {'varexplained', 'mask', 'recomp_beta_scramble.mgz', 'x_scramble.mgz', 'y_scramble.mgz', 'sigma_scramble.mgz'};
elseif (~opt.useBensonMaps && opt.useSmoothedData)
    prfParams = {'varexplained', 'mask', 'recomp_beta', 'x_smoothed', 'y_smoothed', 'sigma_smoothed'};
elseif opt.useBensonMaps
    prfParams = {'mask', 'beta', 'x', 'y', 'sigma'};
else
    prfParams = {'varexplained', 'mask', 'x', 'y', 'sigma', 'beta'};
end


% Load prf parameters
prf = loadpRFsfromSurface(prfParams, prfSurfPath, opt);

% If perturb original pRFs, check dimensions with loaded pRF data
if strcmp(opt.perturbOrigPRFs, 'position')
    assert(size(prf.x_vary,2)==length(opt.varyPosition))
    nIter = length(opt.varyPosition);
elseif strcmp(opt.perturbOrigPRFs, 'size')
    assert(size(prf.sigma_vary,2)==length(opt.varySize))
    nIter = length(opt.varySize);
elseif strcmp(opt.perturbOrigPRFs, 'scramble')
    assert(size(prf.sigma_scramble,2)==opt.nScrambles)
    nIter = opt.nScrambles;
elseif ~opt.perturbOrigPRFs
    nIter = 1;
end

% Keep a copy of all parameters
prfAll = prf; 

% Remove blink periods from stim,
conditions = stim.conditions;
blinkIm    = conditions.triggers.stimConditions(1:size(stim.im,2))==20;
stim.im(:,blinkIm) = NaN;

% Define time/epoch dimension for allocating space and plotting
t = (0:size(stim.im,2)-1) .* diff(opt.epochStartEnd);

% Allocate space
predSurfResponse = NaN(length(t),size(prf.roimask,1),nIter);

% loop over dimensions, if necessary
for ii = 1:nIter
    
    % Select new prf parameters, if they vary in size or position
    if strcmp(opt.perturbOrigPRFs,'position')
        prf.x_vary = prfAll.x_vary(:,ii);
        prf.y_vary = prfAll.y_vary(:,ii);
    elseif strcmp(opt.perturbOrigPRFs,'size')
        prf.sigma_vary = prfAll.sigma_vary(:,ii);
    elseif strcmp(opt.perturbOrigPRFs,'scramble')
        prf.x_scramble = prfAll.x_scramble(:,ii);
        prf.y_scramble = prfAll.y_scramble(:,ii);
        prf.sigma_scramble = prfAll.sigma_scramble(:,ii);
        prf.recomp_beta_scramble = prfAll.recomp_beta_scramble(:,ii); 
    end
    
    % Get predicted response from prf data
    predSurfResponse(:,:,ii) = mprf_MEGPredictionFromSurface(prf, stim); 
end

% Plot predicted response BS surface
if opt.verbose
    if opt.saveFig && ~exist(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'predSurfResponse'),'dir')
         mkdir(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'predSurfResponse')); 
    end
    
    figure(1), set(gcf, 'Position', [652   784   908   554], 'Color', 'w'); 
    for ii = 1:nIter
        clf;
        plot(t,predSurfResponse(:,:,ii));
        title(sprintf('Predicted response to retinotopy stim using pRFs from cortex, %d', ii))
        xlabel('Time (s)'); ylabel('Predicted vertex response (a.u.)');
        set(gca, 'FontSize', 14, 'TickDir','out'); box off
    
        if opt.saveFig
            print(fullfile(dirPth.model.saveFigPth, opt.subfolder, ...
             'predSurfResponse', sprintf('predPRFResponseFromSurface%s_%d', opt.fNamePostFix, ii)), '-dpng')
        end
    end
end

% Remove last dimension out, if not used
predSurfResponse = squeeze(predSurfResponse);

% Save predicted response to MEG stimuli for every vertex
if opt.doSaveData
    if ~exist(fullfile(dirPth.model.saveDataPth,opt.subfolder, 'pred_resp'), 'dir')
        mkdir(fullfile(dirPth.model.saveDataPth, opt.subfolder, 'pred_resp')); end
    save(fullfile(dirPth.model.saveDataPth,opt.subfolder,'pred_resp','predSurfResponseFromPRFs'), 'predSurfResponse', '-v7.3');
end




return