function [prf, prfSurfPath] = mprf_varyPRFPositionOnSurface(prfSurfPath, opt)

if opt.verbose; fprintf('(%s): Vary pRF position...\n', mfilename); end

% Get range to vary prf position
if (~isfield(opt,'varyPosition') || isempty(opt.varyPosition))
    error('(%s): Range to vary prf position is not defined in opt', mfilename)
else
    range = opt.varyPosition;
end

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
    surfmatfname = fullfile(prfSurfPath,'preturbed_prf_params.mat');
    save(surfmatfname, 'prf');
end

return