function [prfData, prfDataPath] = loadBensonRetinotopyMaps(subjID, dirPth, opt)

% Check if freesurfer matlab toolbox is added
if ~exist('MRIread', 'file')
    fshome = getenv('FREESURFER_HOME');
    if isempty(fshome)
        setenv('FREESURFER_HOME', '/Applications/freesurfer');
        fshome = getenv('FREESURFER_HOME');
    end
    addpath(genpath(fullfile(fshome,'matlab')));
    addpath(genpath(fullfile(fshome,'fsfast','toolbox')));
end

[bsDB,bsProtocol] = fileparts(dirPth.bsPth);
bsAnatDir = dirPth.bs.anatPth;
fsdir     = dirPth.fsPth;

% Check if we use full size or downsampled surfaces and gain matrix
if opt.fullSizeGainMtx
    highres = true;
    subFolder = 'highres';
else
    highres = false;
    subFolder = 'lowres';
end

% create folder to save figures
saveDir  = fullfile(dirPth.fmri.saveDataPth, 'BensonAtlasPred', subFolder);
if ~exist(saveDir, 'dir'); mkdir(saveDir); end


%% 1.Create combined hemi template that matches vertices in BS Gain matrix

% Note: time consuming and computationally heavy if you want a high resolution template!
if (opt.fullSizeGainMtx) && (~exist(fullfile(bsAnatDir, 'highres', 'benson14areas_overlay.mat'), 'file'))
    interp_retinotopy(bsDB, fsdir, subjID, subjID, bsProtocol, highres)
end


if (~opt.fullSizeGainMtx) && (~exist(fullfile(bsAnatDir, 'lowres', 'benson14areas_overlay.mat'), 'file'))
    interp_retinotopy(bsDB, fsdir, subjID, subjID, bsProtocol)
end


%% 2. Load retinotopy templates

prfData = struct();

% Load V1-3 template with unitized phases in downsampled brainstorm format (computed by interp_retinotopy.m)
areas    = load(fullfile(bsAnatDir, subFolder,'benson14areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex. CHECK: Positive values represent LH (?) Negative values RH (?)
eccen    = load(fullfile(bsAnatDir, subFolder,'benson14eccen_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred eccentricity in degrees, zeros refer to outside of visual cortex
polarang = load(fullfile(bsAnatDir, subFolder,'benson14angle_overlay.mat')); % [1xNsources] Values represents vertex preferred polar angle in degrees (0 is the upper vertical meridian), zeros refer to outside of visual cortex
sigma    = load(fullfile(bsAnatDir, subFolder,'benson14sigma_overlay.mat')); % [1xNsources] Values represents vertex preferred polar angle in degrees (0 is the upper vertical meridian), zeros refer to outside of visual cortex

% get areas
areas.v1 = (areas.sub_bs_areas==1);
areas.v2 = (areas.sub_bs_areas==2);
areas.v3 = (areas.sub_bs_areas==3);
areas.all = areas.sub_bs_areas>0;

%% 3. Get sigma, x,y coords from polar angle and eccen, and beta from sigma

% note, RH hemi should have negative polar angle values when converting to
% x, y coordinates
% theta = pi/180 * (90 - polarang.sub_bs_angle);
eccenThresh = ((eccen.sub_bs_eccen > 0.5) & (eccen.sub_bs_eccen < opt.eccThresh(2)));
eccentricity =  eccen.sub_bs_eccen .* eccenThresh;

x_tmp = eccentricity .* cos(polarang.sub_bs_angle);
y_tmp = eccentricity .* sin(polarang.sub_bs_angle);

prfData.y            = y_tmp;
prfData.x            = x_tmp;
prfData.eccentricity = eccentricity;
prfData.polar_angle  = polarang.sub_bs_angle;
prfData.sigma        = sigma.sub_bs_sigma;
prfData.beta         = 1./sqrt(2*pi*prfData.sigma.^2);
prfData.mask         = areas.all;

% Store path where files are saved
fn = fieldnames(prfData);
for ii = 1:length(fn)
    fname = fullfile(saveDir, ['pial.',fn{ii}]);
    write_curv(fname,prfData.(fn{ii}),1);
end
    
% Rename path for output
prfDataPath = saveDir;

return
