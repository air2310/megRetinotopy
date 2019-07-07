function predResponseAllVertices = mprf_MEGPredictionFromSurface(prfSurfPath, stim, dirPth, opt)

% Function to make predictions for every vertex on given surface
% 
%   [pred_response, model] = mprf_MEGPredictionFromSurface(prfSurfPath, meg.stim, dirPth, opt);
%
% INPUTS:
%   prfSurfPath     : path to prf parameters on surface (string)
%   stim            : meg stimulus struct (should contain tktktktk)
%   dirPth          : paths to files for given subject
%   opt             : struct with boolean flags
%
% OUTPUTS:
%   prf        : predicted response time series for every vertex on brainstorm
%                surface

%% Find prf parameter files

% Define pRF parameters to read from surface
if (~opt.useBensonMaps && opt.useSmoothedData)
    prfParams = {'varexplained', 'mask', 'recomp_beta', 'x_smoothed', 'y_smoothed', 'sigma_smoothed', };
elseif opt.useBensonMaps
    prfParams = {'mask', 'beta', 'x', 'y', 'sigma'};
else
    prfParams = {'varexplained', 'mask', 'x', 'y', 'sigma', 'beta'};
end

% Get meg stim
stim = stim.meg_stim;

% Remove empty files 
d = dir(fullfile(prfSurfPath, '*'));
for ii = 1:length(d)
    if d(ii).bytes<1
        emptyFile(ii) = 1;
    else
        emptyFile(ii) = 0;
    end
end
d(find(emptyFile)) = []; %#ok<FNDSB>

% Display prf parameters
fprintf('(%s): Found the following prf parameters on the surface: \n',  mfilename)
fprintf('\t %s \n', d.name)


%% Get prf parameters
prf = struct();
for idx = 1:length(prfParams)
    
    % 1. Load prf params
    param = dir(fullfile(prfSurfPath, sprintf('*.%s',prfParams{idx})));
    theseData = read_curv(fullfile(prfSurfPath,param.name));
    
    switch prfParams{idx}
        
        case 'varexplained'
            prf.(prfParams{idx}) = theseData;
            
            % 2. Check variance explained by pRF model and make mask if requested
            if any(opt.varExplThresh)    
                prf.vemask = ((prf.varexplained > opt.varExplThresh(1)) & (prf.varexplained < opt.varExplThresh(2)));
            else % 3. If not, make mask with all ones
                prf.vemask = true(size(prf.varexplained));
            end
    
    
        case 'mask' % 4. Make roi mask (original file: NaN = outside mask, 0 = inside mask)
            prf.roimask = ~isnan(theseData);
        
            % Benson maps don't have variance explained map, so we just use the
            % same vertices as the roi mask
            if opt.useBensonMaps; prf.vemask = prf.roimask; end
            
        case {'beta', 'recomp_beta'}
            if any(opt.betaPrctileThresh)
                thresh = prctile(theseData, opt.betaPrctileThresh);
                betamask = ((theseData > thresh(1)) & (theseData < thresh(2)));
                theseData(~betamask) = NaN;
                prf.(prfParams{idx}) = theseData(prf.vemask & prf.roimask);
            else
                prf.(prfParams{idx}) = theseData(prf.vemask & prf.roimask);
            end
    
        case {'x_smoothed', 'x', 'y_smoothed', 'y', 'sigma_smoothed', 'sigma'}
            prf.(prfParams{idx}) = theseData(prf.vemask & prf.roimask);
    end
            
end

%% Get RFs from prf surface parameters (stim locations x vertices)
fn = fieldnames(prf);
x0    = prf.(fn{cellfind(regexp(fn, '\<x'))});
y0    = prf.(fn{cellfind(regexp(fn, '\<y'))});
sigma = prf.(fn{cellfind(regexp(fn, '\<sigma'))});
beta  = prf.(fn{cellfind(regexp(fn, 'beta'))});

% Build RFs that fall within the stimulus aperture
RF = rfGaussian2d(stim.X,stim.Y,sigma,sigma,false, x0, y0);

% Get predicted response from RFs given stimulus (epochs x vertices)
predResponse = bsxfun(@times, stim.im' * RF, beta');

predResponseAllVertices = NaN(size(predResponse,1), size(prf.vemask,1));
predResponseAllVertices(:, prf.vemask & prf.roimask) = predResponse;

% Plot predicted response BS surface
if opt.verbose
    t = (0:size(stim.im,2)-1) .* diff(opt.epochStartEnd);
    figure, set(gcf, 'Position', [652   784   908   554], 'Color', 'w');
    plot(t,predResponseAllVertices);
    title('Predicted response to stimulus from all vertices, using pRF parameters from surface')
    xlabel('Time (s)'); ylabel('Predicted vertex response (a.u.)');
    set(gca, 'FontSize', 14, 'TickDir','out'); box off
     if opt.saveFig; 
         if ~exist(dirPth.model.saveFigPth, 'dir'); mkdir(dirPth.model.saveFigPth); end
         print(fullfile(dirPth.model.saveFigPth, ...
            sprintf('predPRFResponseFromSurface_benson%d_highres%d_smoothed%d', ...
            opt.useBensonMaps, opt.fullSizeGainMtx, opt.useSmoothedData)), '-dpng')
    end
end

% save predicted response to MEG stimuli for every vertex
if opt.doSaveData
    if ~exist(fullfile(prfSurfPath,'pred_resp'), 'dir'); mkdir(fullfile(prfSurfPath,'pred_resp')); end;
    save(fullfile(prfSurfPath,'pred_resp','predResponseAllVertices'));
end
