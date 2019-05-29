%% s_BensonAtlasPRFPredictions

if ~exist('MRIread') 
    addpath(genpath('/Applications/freesurfer/matlab'));
    addpath(genpath('/Applications/freesurfer/fsfast/toolbox'))
end

% Define params
subject = 'wlsubj040';
exp     = 'MEG_Retinotopy';

% Freesurfer, brainstorm database and data/anatomy directories
fsdir    = '/Volumes/server/Freesurfer_subjects';
bsDB     = '/Volumes/server/Projects/MEG/brainstorm_db';
megRet   = '/Volumes/server/Projects/MEG/Retinotopy/Subject_sessions/';

d = dir(fullfile(bsDB, exp, 'data', subject));
datadir = fullfile(bsDB, exp, 'data', subject, d(end).name);
anatdir = fullfile(bsDB, exp, 'anat', subject);

%% 1.1 Create downsampled amplitude template
interp_retinotopy(bsDB, fsdir, subject, subject, exp)

%% 1. Load Gain matrix
keep_sensors    = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

G_constrained = getGainMatrix(datadir, keep_sensors);

%% 2. Load retinotopy templates (downsampled)

% Load V1-3 template with unitized phases in downsampled brainstorm format (computed by interp_retinotopy.m)
areas    = load(fullfile(anatdir, 'benson14_varea.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex. CHECK: Positive values represent LH (?) Negative values RH (?)
eccen    = load(fullfile(anatdir, 'benson14_eccen.mat')); % [1xNsources] Nonzero value represents vertex preferred eccentricity in degrees, zeros refer to outside of visual cortex
polarang = load(fullfile(anatdir, 'benson14_angle.mat')); % [1xNsources] Values represents vertex preferred polar angle in degrees (0 is the upper vertical meridian), zeros refer to outside of visual cortex
sigma    = load(fullfile(anatdir, 'benson14_sigma.mat')); % [1xNsources] Values represents vertex preferred polar angle in degrees (0 is the upper vertical meridian), zeros refer to outside of visual cortex

stim  = load(fullfile(megRet, subject, 'stimuli', 'meg', 'imported_stimulus', 'meg_stimulus.mat'));

% RF matrix:
    RF = rfGaussian2d(stim.X,stim.Y,cur_sigma,cur_sigma,false, cur_x0,cur_y0);
    
    % Predicted BOLD responses:
    tmp = stim.im_conv' * RF;