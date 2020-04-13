%% s_initSession_megRet

% This is a script to define the paths of the files we need to run a model.
% 

% Define subject and its paths
subject            = 'wlsubj058';
dataDir            = '/Volumes/server/Projects/MEG/Retinotopy/Data/';
freeSurferDir      = '/Volumes/server/Freesurfer_subjects/';
brainstormDBDir    = '/Volumes/server/Projects/MEG/brainstorm_db/';
vistaSessionDir    = fullfile(dataDir, 'fMRI', subject, 'vistaSession');
megDataDir         = fullfile(dataDir, 'MEG', subject);
outPutDir          = fullfile('/Volumes/server/Projects/MEG/Retinotopy/', 'Subject_sessions', subject);

% Predefine struct
s = struct();

% Add FS and BS subject
s.fsSubject          = subject;
s.bsSubject          = subject;

% Find MEG time series
s.MEGtimeseries = fullfile(megDataDir, 'processed', 'epoched_data_hp_preproc_denoised.mat');

% Add other general paths
s.vistaSession.pth = vistaSessionDir;
s.freeSurferDir    = freeSurferDir;
s.megData.pth      = megDataDir;
s.outPut.pth       = outPutDir;

% Find MRI PRF params
d = dir(fullfile(vistaSessionDir, 'Gray', '*','*fFit*'));
s.PRFParams.pth    = fullfile(d.folder, d.name);

% Find MEG & MRI stimulus (and check if they are the same)
s.MRIStimIm.pth     = fullfile(vistaSessionDir, 'Stimuli', 'scan_images.mat');
s.MRIStimParams.pth = fullfile(vistaSessionDir, 'Stimuli', 'scan_params.mat');
s.MEGStim.pth       = fullfile(megDataDir, 'stimFiles', 'MEG_retinotopy_stimulus_run_1.mat');
s.MEGStimGrid.pth   = fullfile(megDataDir, 'stimFiles', 'MEG_grid.mat');

% Find BS surfaces
s.BS.surface.pth = fullfile(brainstormDBDir, 'MEG_Retinopy', 'anat', subject);

% Find Brainstorm headmodel and copy to output dir
d = dir(fullfile(brainstormDBDir, 'MEG_Retinotopy', 'data', subject, '*', 'headmodel*.mat'));
s.BS.gainMatrix.pth = fullfile(d.folder,d.name);

%% Find ROIS in Freesurfer directory
s.ROIs.pth = fullfile(freeSurferDir, subject, 'surf' );

%% Copy gainmatrix to subject folder
status = copyfile(s.BS.gainMatrix.pth, fullfile(outPutDir, d.name));

%% Load stimulus (To Do: check where coordinates fall wrt image)
s.stim = loadStim(s); 

%% Smooth prf parameters
% megRet_smoothPRFParamsCorrected(s)
megRet_smoothPRFParams(s)

% Go from FS to BS
roiType = 'allRoisWangAtlas'; % What rois?
megRet_FS2BS(s, roiType)





