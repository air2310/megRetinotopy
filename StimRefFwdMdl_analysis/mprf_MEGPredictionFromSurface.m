function predResponseAllVertices = mprf_MEGPredictionFromSurface(prfSurfPath, stim, subjID, dirPth, opt)

% Function to make predictions for every vertex on given surface
% 
%   [pred_response, model] = mprf_MEGPredictionFromSurface(prfSurfPath, meg.stim, subjID, dirPth, opt);
%
% INPUTS:
%   prfSurfPath     : path to prf parameters on surface (string)
%   stim            : meg stimulus struct (should contain tktktktk)
%   subjID          : subject name (string)
%   dirPth          : paths to files for given subject
%   opt             : struct with boolean flags
%
% OUTPUTS:
%   prf        : predicted response time series for every vertex on brainstorm
%                surface

%% Find prf parameter files

% Define pRF parameters to read from surface
if opt.useSmoothedData
    prfParams = {'varexplained', 'mask', 'x_smoothed', 'y_smoothed', 'sigma_smoothed', 'recomp_beta'};
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
d(find(emptyFile)) = [];

% Display prf parameters
sprintf('Found surface prf parameters: %s \n', d.name)


%% Get prf parameters
prf = struct();
for idx = 1:length(prfParams)
    
    % 1. Load prf params
    param = dir(fullfile(prfSurfPath, sprintf('*.%s',prfParams{idx})));
    theseData = read_curv(fullfile(prfSurfPath,param.name));
    
    % 2. Check variance explained by pRF model and make mask if requested
    if ((strcmp(prfParams{idx},'varexplained')) && any(opt.varExplThresh))
        prf.(prfParams{idx}) = theseData;
        prf.vemask = ((prf.varexplained > opt.varExplThresh(1)) & (prf.varexplained < opt.varExplThresh(2)));
    
    % 3. If not, make mask with all ones
    elseif strcmp(prfParams{idx},'varexplained')
        prf.(prfParams{idx}) = theseData;
        prf.vemask = true(size(prf.varexplained));
    
    % 4. Make roi mask
    elseif strcmp(prfParams{idx},'mask') % (original file: NaN = outside mask, 0 = inside mask)
        prf.roimask = ~isnan(theseData);   
    
    % 5. Mask data
    else  
        prf.(prfParams{idx}) = theseData(prf.vemask & prf.roimask);        
    end

end

%% Get RFs from prf surface parameters (stim locations x vertices)
x0    = prf.(prfParams{3});
y0    = prf.(prfParams{4});
sigma = prf.(prfParams{5});
beta  = prf.(prfParams{6});

% Build RFs that fall within the stimulus aperture
RF = rfGaussian2d(stim.X,stim.Y,sigma,sigma,false, x0, y0);

% Get predicted response from RFs given stimulus (epochs x vertices)
predResponse = bsxfun(@times, stim.im' * RF, beta');

predResponseAllVertices = NaN(size(predResponse,1), size(prf.vemask,1));
predResponseAllVertices(:, prf.vemask & prf.roimask) = predResponse;

% Plot predicted response BS surface
if opt.verbose
    figure, plot(1:140,predResponseAllVertices)
end

% save predicted response to MEG stimuli for every vertex
if opt.doSaveData
    if ~exist(fullfile(prfSurfPath,'pred_resp')); mkdir(fullfile(prfSurfPath,'pred_resp')); end;
    save(fullfile(prfSurfPath,'pred_resp'),'predResponseAllVertices');
end
