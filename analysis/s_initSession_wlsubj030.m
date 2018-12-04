%% s_initSession_wlsubj030

% This is a script to define the paths of the files we need to run a model.
% 

% Define subject and its paths
subject            = 'wlsubj030';
dataDir            = '/Volumes/server/Projects/MEG/Retinotopy/Data/';
freeSurferDir      = '/Volumes/server/Freesurfer_subjects/';
brainstormDBDir    = '/Volumes/server/Projects/MEG/brainstorm_db/';
vistaSessionDir    = fullfile(dataDir, 'fMRI', subject, 'VistaSession');
megDataDir         = fullfile(dataDir, 'MEG', subject);
outPutDir          = fullfile(mprfRootPath, 'data', 'subjectSession', subject);

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
s.MEGStim.pth       = fullfile(megDataDir,'stimFiles', 'MEG_retinotopy_stimulus_run_1.mat');
s.MEGStimGrid.pth   = fullfile(megDataDir,'stimFiles', 'MEG_grid.mat');
s.BS.surface.pth = fullfile(brainstormDBDir, 'MEG_Retinopy', 'anat', subject);

%% Find Brainstorm headmodel and copy to output dir
d = dir(fullfile(brainstormDBDir, 'MEG_Retinopy', 'data', subject, '*', 'headmodel*.mat'));
s.BS.gainMatrix.pth = fullfile(d.folder,d.name);
status = copyfile(s.BS.gainMatrix.pth, fullfile(outPutDir, d.name));

%% Load stimulus (To Do: check where coordinates fall wrt image)
s.stim = loadStim(s); 

%% Smooth prf parameters
megRet_smoothPRFParams(s)

%% Transform FS 2 BS 

% Find ROIS in Freesurfer directory
s.ROIs.pth = fullfile(freeSurferDir, subject, 'surf', 'WangIndividualROIs');

% Go from FS to BS
roiType = 'allRoisWangAtlas'; % What rois?
megRet_FS2BS(s, roiType)



