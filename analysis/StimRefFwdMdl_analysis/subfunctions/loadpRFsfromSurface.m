function prf = loadpRFsfromSurface(prfParams, prfSurfPath, opt)
% Function to vary pRF centers on the cortical surface (rotation in polar
% angle)
%       prf = loadpRFsfromSurface(prfParams, prfSurfPath, opt)
% 
% INPUTS:
%   prfSurfPath     : path to surface files containing prf parameters (str)
%   dirPth          : paths to files for given subject (struct)
%   opt             : struct with boolean flags. Should contain the field 
%                     'perturbOrigPRFs' with one of following definitions:
%                     'position', 'size', 'scramble', or false to exit
%
% OUTPUTS:
%   prf             : struct with prf data on surface found in prfSurfPath 
%
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019



% Check files in folder and remove empty files
d = dir(fullfile(prfSurfPath, '*'));

% Display prf parameters
if opt.verbose
    fprintf('(%s): Found the following prf parameters on the surface: \n',  mfilename)
    fprintf('\t %s \n', d.name)
end

% Create a new struct to add prf data
prf = struct();

% Find X and Y prf parameters
filenameX = ['pial.' prfParams{find(~cellfun(@isempty,regexp(prfParams, '\<x')))}]; %
filenameY = ['pial.' prfParams{find(~cellfun(@isempty,regexp(prfParams, '\<y')))}];

if regexp(filenameX, '.mgz', 'ONCE')
    if opt.mri.useSmoothedData
        x0_tmp = read_curv(fullfile(prfSurfPath,'pial.x_smoothed'));
        y0_tmp = read_curv(fullfile(prfSurfPath,'pial.y_smoothed'));
    else
        x0_tmp = read_curv(fullfile(prfSurfPath,'pial.x'));
        y0_tmp = read_curv(fullfile(prfSurfPath,'pial.y'));
    end
    [~,eccen] = cart2pol(x0_tmp,y0_tmp);
    prf.eccenmask = (eccen<=opt.mri.eccThresh(2));
else
    x0_tmp = read_curv(fullfile(prfSurfPath,filenameX));
    y0_tmp = read_curv(fullfile(prfSurfPath,filenameY));
    [~,eccen] = cart2pol(x0_tmp,y0_tmp);
    prf.eccenmask = (eccen<=opt.mri.eccThresh(2));
end


for idx = 1:length(prfParams)
    
    % Load prf params
    param     = dir(fullfile(prfSurfPath, sprintf('pial.%s',prfParams{idx})));
    
    if regexp(param.name, '.mgz', 'ONCE')
        tmp = MRIread(fullfile(prfSurfPath,param.name));
        theseData = tmp.vol;
    else
        theseData = read_curv(fullfile(prfSurfPath,param.name));
    end

    switch prfParams{idx}
        
        case 'varexplained'
            
            prf.(prfParams{idx}) = theseData;
            
            % Check variance explained by pRF model and make mask if requested
            if any(opt.mri.varExplThresh)
                prf.vemask = ((prf.varexplained >= opt.mri.varExplThresh(1)) & (prf.varexplained < opt.mri.varExplThresh(2)));
            else % If not, make mask with all ones
                prf.vemask = true(size(prf.varexplained));
            end
            
        case {'mask', 'V123mask'} % Make roi mask (original file: 0 = outside mask, 1 = inside mask)
            prf.roimask = (theseData>0);
            
            % Benson maps don't have variance explained map, so we just use the
            % same vertices as the roi mask
            if opt.mri.useBensonMaps; prf.vemask = prf.roimask; end
            
        case {'beta', 'recomp_beta'}
            if any(opt.mri.betaPrctileThresh)
                thresh = prctile(theseData, opt.mri.betaPrctileThresh);
                betamask = ((theseData > thresh(1)) & (theseData < thresh(2)));
                theseData(~betamask) = NaN;
                prf.(prfParams{idx}) = theseData(prf.vemask & prf.roimask & prf.eccenmask);
            else
                prf.(prfParams{idx}) = theseData(prf.vemask & prf.roimask & prf.eccenmask);
            end
            
        case {'x_smoothed', 'x', 'y_smoothed', 'y', 'sigma_smoothed', 'sigma'}
            prf.(prfParams{idx}) = theseData(prf.vemask & prf.roimask & prf.eccenmask);
            
        case {'x_vary.mgz', 'y_vary.mgz', 'sigma_vary.mgz', ...
                'x_smoothed_vary.mgz', 'y_smoothed_vary.mgz', 'sigma_smoothed_vary.mgz', ...
                'x_smoothed_scramble.mgz', 'y_smoothed_scramble.mgz', 'sigma_smoothed_scramble.mgz', ...
                'x_scramble.mgz', 'y_scramble.mgz', 'sigma_scramble.mgz', 'recomp_beta_scramble.mgz'}
            fn = strsplit(prfParams{idx}, '.');
            prf.(fn{1}) = theseData; % no mask needed, since they are created from masked data
    end
end

return
