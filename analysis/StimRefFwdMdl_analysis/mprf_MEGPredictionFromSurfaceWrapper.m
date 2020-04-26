function predSurfResponse = mprf_MEGPredictionFromSurfaceWrapper(prfSurfPath, stim, dirPth, opt)
% Wrapper function for mprf_MEGPredictionFromSurface. This wrapper was
% created to keep the function to make predictions from prf data simple,
% while still allowing multiple analysis options.
%
% INPUTS:
%   prfSurfPath     : directory where prf data live (string)
%   stim            : meg stimulus struct (should contain fieldname 'im',
%                       with windowed stimulus images for every epoch)
%   dirPth          : paths to files for given subject (struct)
%   opt             : struct with boolean flags
%
% OUTPUTS:
%   predSurfResponse : predicted response time series (vertex x timepoints [x optional variations])
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019

fprintf('(%s): Make predictions on cortical surface.', mfilename)

% Create folders if figures are saved
if opt.saveFig
    if ~exist(fullfile(dirPth.model.saveFigPth, opt.subfolder,'predSurfResponse'),'dir')
        mkdir(fullfile(dirPth.model.saveFigPth, opt.subfolder,'predSurfResponse'));
    end
end

% Define pRF parameters to read from surface
prfParams = getpRFParamNames(opt);

% Load prf parameters
prf = loadpRFsfromSurface(prfParams, prfSurfPath, opt);

% Check dimensions with loaded pRF data, and set the number of iterations
iter  = checkNumberOfIterations(prf, opt, 'prfSurf');
nIter = length(iter);

% Keep a copy of all parameters
prfAll = prf;

% Remove blink periods from stim,
conditions = stim.conditions;
blinkIm    = conditions.stimConditions(1:size(stim.im,2))==20;
stim.im(:,blinkIm) = NaN;

% Define time/epoch dimension for allocating space and plotting
t = (0:size(stim.im,2)-1) .* diff(opt.meg.epochStartEnd);

% Allocate space
predSurfResponse = NaN(length(t),size(prf.varexplained,1),nIter);


%% loop over dimensions, if necessary
for ii = 1:nIter
    
    fprintf('.')
    % Select new prf parameters, if they vary in size or position
    if strcmp(opt.vary.perturbOrigPRFs,'position')
        if opt.mri.useSmoothedData
            prf.x_smoothed_vary     = prfAll.x_smoothed_vary(:,ii);
            prf.y_smoothed_vary     = prfAll.y_smoothed_vary(:,ii);
        else
            prf.x_vary              = prfAll.x_vary(:,ii);
            prf.y_vary              = prfAll.y_vary(:,ii);
        end
    elseif strcmp(opt.vary.perturbOrigPRFs,'size')
        if opt.mri.useSmoothedData
            prf.sigma_smoothed_vary = prfAll.sigma_smoothed_vary(:,ii);
        else
            prf.sigma_vary          = prfAll.sigma_vary(:,ii);
        end
    elseif strcmp(opt.vary.perturbOrigPRFs,'scramble')
        if opt.mri.useSmoothedData
            prf.x_smoothed_scramble     = prfAll.x_smoothed_scramble(:,ii);
            prf.y_smoothed_scramble     = prfAll.y_smoothed_scramble(:,ii);
            prf.sigma_smoothed_scramble = prfAll.sigma_smoothed_scramble(:,ii);
            prf.recomp_beta_scramble    = prfAll.recomp_beta_scramble(:,ii);
        else
            prf.x_scramble          = prfAll.x_scramble(:,ii);
            prf.y_scramble          = prfAll.y_scramble(:,ii);
            prf.sigma_scramble      = prfAll.sigma_scramble(:,ii);
            prf.beta_scramble       = prfAll.beta_scramble(:,ii);
        end
    end
    
    % MAIN FUNCTION: Get predicted response from prf data
    predSurfResponse(:,:,ii) = mprf_MEGPredictionFromSurface(prf, stim);
end

% Check original location if perturbing prfs, we need this iteration to 
% base off the variance mask
if strcmp(opt.vary.perturbOrigPRFs,'position')
    origPRFidx = find(opt.vary.position==0);
elseif strcmp(opt.vary.perturbOrigPRFs,'size')
    origPRFidx = find(opt.vary.size==1);
elseif strcmp(opt.vary.perturbOrigPRFs,'scramble')
    origPRFidx = 1;
else
    origPRFidx = 1;
end

% Check for outliers in variance
if opt.mri.predSurfMaxThresh(2)
    origPRFPredSurfResponse = predSurfResponse(:,:,origPRFidx);
    maxSurfResp   = max(origPRFPredSurfResponse,[],'omitnan');
    grand_median  = nanmedian(maxSurfResp(maxSurfResp>0));
    exclOrigVerts = (maxSurfResp > opt.mri.predSurfMaxThresh(2)*grand_median);
    predSurfResponse(:,exclOrigVerts,:) = NaN;
end

fprintf('.Done!\n');

%% Plot predicted response BS surface
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
if opt.saveData
    if ~exist(fullfile(dirPth.model.saveDataPth,opt.subfolder, 'pred_resp'), 'dir')
        mkdir(fullfile(dirPth.model.saveDataPth, opt.subfolder, 'pred_resp')); end
    save(fullfile(dirPth.model.saveDataPth,opt.subfolder,'pred_resp','predSurfResponseFromPRFs'), 'predSurfResponse', '-v7.3');
end

% Create video of responses with stimulus on the side
if (opt.verbose && any(~opt.vary.perturbOrigPRFs))
    %    makeSurfResponseMovie(predSurfResponse, stim, dirPth, opt);
end


return
