%% s_initSession_wlsubj030


% Define subject and its paths
subject            = 'wlsubj030';
dataDir            = '/Volumes/server/Projects/MEG/Retinotopy/Data/';
vistaSessionDir    = fullfile(dataDir, 'fMRI', subject, 'VistaSession');
megPreprocessedDir = fullfile(dataDir, 'MEG', subject, 'processed');
outPutDir          = fullfile('/Volumes/server/Projects/MEG/Subject_sessions', subject);

s = struct();

% Find MEG time series
s.MEGtimeseries = fullfile(megPreprocessedDir, 'epoched_data_hp_preproc_denoised.mat');

% Find MRI PRF params
s.PRFParams = fullfile(vistaSessionDir, 'Gray', 'rm_);


% Find stimulus

% Find ROIS in Freesurfer directory

% Find Brainstorm headmodel

% Export FS surfaces to Brainstorm 

