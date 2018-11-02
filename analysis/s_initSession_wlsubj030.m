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

% Add other general paths
s.vistaSession.pth = vistaSessionDir;
s.megData.pth      = megDataDir;
s.outPut.pth       = outPutDir;

% Find MRI PRF params
d = dir(fullfile(vistaSessionDir, 'Gray', '*','*fFit*'));
s.PRFParams.pth    = fullfile(d.folder, d.name);

% Find MEG & MRI stimulus (and check if they are the same)
s.MRIStimIm.pth     = fullfile(vistaSessionDir, 'Stimuli', 'scan_images.mat');
s.MRIStimParams.pth = fullfile(vistaSessionDir, 'Stimuli', 'scan_params.mat');
s.MEGStim.pth       = fullfile(megDataDir, 'raw', 'R0942_MegRet_9.11.18','stimFiles', 'MEG_retinotopy_stimulus_run_1.mat');
s.MEGStimGrid.pth   = fullfile(megDataDir, 'raw', 'R0942_MegRet_9.11.18','stimFiles', 'MEG_grid.mat');

% Load stimulus (To Do: check where coordinates fall wrt image)
s.stim = loadStim(s); 

% Smooth prf parameters
megRet_smoothPRFParams(s)

% Find ROIS in Freesurfer directory

% Find Brainstorm headmodel

% Export FS surfaces to Brainstorm 

