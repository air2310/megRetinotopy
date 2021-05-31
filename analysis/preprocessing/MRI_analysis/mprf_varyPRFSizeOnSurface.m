function prf = mprf_varyPRFSizeOnSurface(prfSurfPath, opt)
% Function to perturb original pRF parameters, by scaling the original pRF 
% size estimated with fMRI on the cortical surface.
%
%       prf = mprf_varyPRFSizeOnSurface(prfSurfPath, opt)
% 
% INPUTS:
%   prfSurfPath     : path to surface files containing prf parameters (string)
%   opt             : struct with boolean flags. Should contain the field 
%                     'perturbOrigPRFs' defined as 'size' and have
%                     a field called 'varySize' with a vector of scale 
%                     factors to use.
%
% OUTPUT:
%   prf             : struct with prf data, separate for every parameter
%
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019

% Get range to vary prf position
if (~isfield(opt.vary,'size') || isempty(opt.vary.size))
    error('(%s): Range to vary prf size is not defined in opt', mfilename)
else
    range = opt.vary.size;
end

if opt.verbose; fprintf('(%s): Scale pRF sizes with factors: %s\n', mfilename, sprintf('%1.1f ', range)); end


% Load prf parameters on surface
if opt.mri.useBensonMaps
    prfParams = {'mask', 'beta','x','y','sigma'};
elseif opt.mri.useHCPAveMaps || opt.mri.useNYU3TAveMaps
    prfParams = {'varexplained', 'mask','beta','x','y','sigma'};
elseif opt.mri.useSmoothedData
    prfParams = {'varexplained', 'mask','recomp_beta','x_smoothed','y_smoothed','sigma_smoothed'};
else
    prfParams = {'varexplained', 'mask','beta','x','y','sigma'};
end

prf = loadpRFsfromSurface(prfParams, prfSurfPath, opt);

% Get fieldnames to add vary sigma data
fn = fieldnames(prf);
sigma    = prf.(fn{cellfind(regexp(fn, '\<sigma'))});
fn_varySigma = strcat(fn{cellfind(regexp(fn, '\<sigma'))},'_vary');

% Multiply original sigma with range of scaling factors
prf.(fn_varySigma) = sigma * range;


% Create file names for x,y vary position data and save in same folder
if opt.saveData
    
    surfdata = prf.(fn_varySigma);
    surfdatamgzfile = sprintf('%s/%s.%s.mgz', prfSurfPath, 'pial', fn_varySigma);
    % Note save as mgz volume, since write_curv only allows vector
    MRIwrite(struct('vol', surfdata), surfdatamgzfile);

    % Save smoothed prf parameters
    surfmatfname = fullfile(prfSurfPath,'perturbed_prf_params_vary_size.mat');
    save(surfmatfname, 'prf');
end

return