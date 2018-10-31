%% s_initSession_wlsubj030

% This is a script to define the paths of the files we need to run a model.
% 

% Define subject and its paths
subject            = 'wlsubj030';
dataDir            = '/Volumes/server/Projects/MEG/Retinotopy/Data/';
vistaSessionDir    = fullfile(dataDir, 'fMRI', subject, 'VistaSession');
megDataDir         = fullfile(dataDir, 'MEG', subject);
outPutDir          = fullfile(mprfRootPath, 'data', subject);

% Predefine struct
s = struct();

% Find MEG time series
s.MEGtimeseries = fullfile(megDataDir, 'processed', 'epoched_data_hp_preproc_denoised.mat');

% Find MRI PRF params
d = dir(fullfile(vistaSessionDir, 'Gray', '*','*fFit*'));
s.PRFParams.pth = fullfile(d.folder, d.name);

% Find MEG & MRI stimulus (and check if they are the same)
s.MRIStimIm.pth     = fullfile(vistaSessionDir, 'Stimuli', 'scan_images.mat');
s.MRIStimParams.pth = fullfile(vistaSessionDir, 'Stimuli', 'scan_params.mat');
s.MEGStim.pth       = fullfile(megDataDir, 'raw', 'R0942_MegRet_9.11.18','stimFiles', 'MEG_retinotopy_stimulus_run_1.mat');
s.MEGStimGrid.pth   = fullfile(megDataDir, 'raw', 'R0942_MegRet_9.11.18','stimFiles', 'MEG_grid.mat');

% RM model, with stimulus
prfParams = load(s.PRFParams.pth);

megStimulus = load(s.MEGStim.pth, 'stimulus');
megStimulusGrid = load(s.MEGStimGrid.pth);

%   Make grid?
sz1 = numel(X);
sz2 = 
X   = repmat(X(:),1,sz2);
Y   = repmat(Y(:),1,sz2);

% define centered around (0,0) coordinate space
x0 = repmat(x0(:)',sz1,1);
y0 = repmat(y0(:)',sz1,1);

% Translate grid so that center is at RF center
X = X - x0;   % positive x0 moves center right
Y = Y - y0;   % positive y0 moves center up
    
    
% check where coordinates fall wrt image?

% Find ROIS in Freesurfer directory

% Find Brainstorm headmodel

% Export FS surfaces to Brainstorm 

