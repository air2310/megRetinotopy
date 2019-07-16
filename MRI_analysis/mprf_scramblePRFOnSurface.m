function prf = mprf_scramblePRFOnSurface(prfSurfPath, opt)
% Function to perturb original pRF parameters, by scrambling the original
% vertices that are linked to a particular pRF estimated with fMRI on the
% cortical surface. We do this by assiging a new vertex location to a pRF
% while keeping the parameters grouped together. The new vertex can be
% anywhere in the ROI mask (so across visual areas, if there are multiple
% visual areas within the ROI mask).
%
%       prf = mprf_scramblePRFOnSurface(prfSurfPath, opt)
%
% INPUTS:
%   prfSurfPath     : path to surface files containing prf parameters (string)
%   opt             : struct with boolean flags. Should contain the field
%                     'perturbOrigPRFs' defined as 'scramble' and a field
%                     called 'nScrambles' with an integer define the number
%                     of scrambling iterations.
%
% OUTPUT:
%   prf             : struct with prf data, separate for every parameter
%
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019

% Get range to vary prf position
if (~isfield(opt,'nScrambles') || isempty(opt.nScrambles))
    error('(%s): Number of scrambling iterations is not defined in opt', mfilename)
end

if opt.verbose; fprintf('(%s): Scramble pRFs  %dx\n', mfilename, opt.nScrambles); end


% Load prf parameters on surface
% if (~opt.useBensonMaps && opt.useSmoothedData)
%     prfParams = {'varexplained', 'recomp_beta', 'mask', 'x_smoothed', 'y_smoothed'};
if opt.useBensonMaps
    prfParams = {'mask', 'beta','x','y','sigma'};
else
    prfParams = {'varexplained', 'mask','recomp_beta','x','y','sigma'};
end

prf = loadpRFsfromSurface(prfParams, prfSurfPath, opt);

% Get fieldnames to add vary sigma data
fn = fieldnames(prf);
x0      = prf.(fn{cellfind(regexp(fn, '\<x'))});
y0      = prf.(fn{cellfind(regexp(fn, '\<y'))});
sigma   = prf.(fn{cellfind(regexp(fn, '\<sigma'))});
beta    = prf.(fn{cellfind(regexp(fn, 'recomp_beta'))});
roimask = prf.(fn{cellfind(regexp(fn, 'roimask'))});
vemask  = prf.(fn{cellfind(regexp(fn, 'vemask'))});


% Update random number generator (change input to integer for reproducibility)
rng('shuffle')

% Get scrambling array
idx = find(roimask&vemask);
scr = randi(length(idx), length(idx), opt.nScrambles);

% Rearrange original prf parameters with new scrambling location
prf.x_scramble            = x0(scr);
prf.y_scramble            = y0(scr);
prf.sigma_scramble        = sigma(scr);
prf.recomp_beta_scramble  = beta(scr);

if opt.doSaveData
    
    % get new filenames
    fn = fieldnames(prf);
    
    for ii = cellfind(regexp(fn, 'scramble'))'
        surfdata = prf.(fn{ii});
        surfdatamgzfile = sprintf('%s/%s.%s.mgz', prfSurfPath, 'pial', fn{ii});
        % Note save as mgz volume, since write_curv only allows vector
        MRIwrite(struct('vol', surfdata), surfdatamgzfile);
        
        % Save smoothed prf parameters
        surfmatfname = fullfile(prfSurfPath,'perturbed_prf_params_scramble.mat');
        save(surfmatfname, 'prf');
    end
end

return