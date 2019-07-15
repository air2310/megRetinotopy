function prf = mprf_varyPRFPositionOnSurface(prfSurfPath, opt)
% Function to perturb original pRF parameters, by rotating the original pRF
% center position estimated with fMRI around the polar angle on the
% cortical surface.
%
%       prf = mprf_varyPRFSizeOnSurface(prfSurfPath, opt)
% 
% INPUTS:
%   prfSurfPath     : path to surface files containing prf parameters (string)
%   opt             : struct with boolean flags. Should contain the field 
%                     'perturbOrigPRFs' defined as 'position' and a field
%                     valled 'varPostion' with a vector that contains the
%                     shift you want the position to rotate along its polar
%                     angle (radians).
%
% OUTPUT:
%   prf             : struct with prf data, separate for every parameter
%
%
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019

% Get range to vary prf position
if (~isfield(opt,'varyPosition') || isempty(opt.varyPosition))
    error('(%s): Range to vary prf position is not defined in opt', mfilename)
else
    range = opt.varyPosition;
end

if opt.verbose; fprintf('(%s): Rotate pRF position along polar angles: %s\n', mfilename, sprintf('%1.1f ', range)); end

% Load prf parameters on surface
% if (~opt.useBensonMaps && opt.useSmoothedData)
%     prfParams = {'varexplained', 'recomp_beta', 'mask', 'x_smoothed', 'y_smoothed'};
if opt.useBensonMaps
    prfParams = {'mask', 'beta','x', 'y'};
else
    prfParams = {'varexplained', 'mask','recomp_beta','x', 'y'};
end

prf = loadpRFsfromSurface(prfParams, prfSurfPath, opt);

% Get fieldnames to locate x and y prf data
fn = fieldnames(prf);
x0    = prf.(fn{cellfind(regexp(fn, '\<x'))});
y0    = prf.(fn{cellfind(regexp(fn, '\<y'))});

fn_varyX = strcat(fn{cellfind(regexp(fn, '\<x'))},'_vary');
fn_varyY = strcat(fn{cellfind(regexp(fn, '\<y'))},'_vary');

prf.(fn_varyX) = (x0*cos(range))-(y0*sin(range));
prf.(fn_varyY) = (x0*sin(range))+(y0*cos(range));


% Create file names for x,y vary position data and save in same folder
if opt.doSaveData
    newFieldNames = {fn_varyX, fn_varyY};
    for ii = 1:length(newFieldNames)
        surfdata = prf.(newFieldNames{ii});
        surfdatamgzfile = sprintf('%s/%s.%s.mgz', prfSurfPath, 'pial', newFieldNames{ii});
        % Note save as mgz volume, since write_curv only allows vector
        MRIwrite(struct('vol', surfdata), surfdatamgzfile);
    end
    % Save smoothed prf parameters
    surfmatfname = fullfile(prfSurfPath,'perturbed_prf_params_vary_pos.mat');
    save(surfmatfname, 'prf');
end

return