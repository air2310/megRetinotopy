function dirPth = loadPaths(subjID)
%
% Function to create a struct with all different directories for the
% different types of files containing data

%% ----- General -----
dirPth  = struct();
dirPth.rootPth = mprf_rootPath;
dirPth.subjID  = subjID;

dirPth.sessionPth     = fullfile(mprf_rootPath, 'data','Retinotopy', 'Subject_sessions', subjID);

%% ------ MEG ------
dirPth.meg.dataPth      = fullfile(mprf_rootPath, 'data','Retinotopy', 'Data', 'MEG'); % should be sym link in git folder
dirPth.meg.saveFigPth   = fullfile(mprf_rootPath, 'data','Retinotopy', 'Quality_check', subjID, 'meg'); % should be sym link in git folder

% Derive other file paths
dirPth.meg.rawSqdPth    = fullfile(dirPth.meg.dataPth, subjID, 'raw');
dirPth.meg.paramFilePth = fullfile(dirPth.meg.dataPth, subjID,  'paramFiles');
dirPth.meg.stimFilePth  = fullfile(dirPth.meg.dataPth, subjID, 'stimFiles');
dirPth.meg.processedDataPth = fullfile(dirPth.meg.dataPth, subjID, 'processed');

dirPth.meg.stimFile     = fullfile(dirPth.meg.stimFilePth, 'MEG_retinotopy_stimulus_run_1.mat');
dirPth.meg.stimGridFile = fullfile(dirPth.meg.stimFilePth, 'MEG_grid.mat');


%% ------ FreeSurfer ------ 
dirPth.fsPth            = fullfile(mprf_rootPath,'data','Freesurfer_subjects'); % should be sym link in git folder

% Derive other file paths
dirPth.fs.segPth        = fullfile(dirPth.fsPth, subjID);
dirPth.fs.surfPth       = fullfile(dirPth.fs.segPth, 'surf');


%% ------ fMRI ------ 
dirPth.fmri.dataPth     = fullfile(mprf_rootPath, 'data','Retinotopy', 'Data', 'fMRI'); % should be sym link in git folder
dirPth.fmri.saveDataPth = fullfile(mprf_rootPath, 'data','Retinotopy', 'Quality_check', subjID, 'fmri'); % should be sym link in git folder

% Derive other file paths
dirPth.fmri.mrvPth       = fullfile(dirPth.fmri.dataPth, subjID, 'vistaSession');
dirPth.fmri.paramFilePth = fullfile(dirPth.fmri.mrvPth, 'Stimuli', 'paramFiles');
dirPth.fmri.stimFilePth  = fullfile(dirPth.fmri.mrvPth, 'Stimuli', 'stimFiles');

dirPth.fmri.saveDataPth_prfMrv = fullfile(dirPth.fmri.saveDataPth, 'prfMrv');
dirPth.fmri.saveDataPth_prfFS = fullfile(dirPth.fmri.saveDataPth, 'prfFS');
dirPth.fmri.saveDataPth_prfBS = fullfile(dirPth.fmri.saveDataPth, 'prfBS');

dirPth.fmri.saveDataPth_roiFS = fullfile(dirPth.fmri.saveDataPth, 'roiFS');
dirPth.fmri.saveDataPth_roiBS = fullfile(dirPth.fmri.saveDataPth, 'roiBS');

d = dir(fullfile(dirPth.fmri.mrvPth, 'Gray', 'Averages','*fFit*'));
dirPth.fmri.vistaGrayFitFile  = fullfile(d.folder, d.name);

%% ------ Brainstorm ------ 
dirPth.bsPth = fullfile(mprf_rootPath, 'data', 'brainstorm_db', 'MEG_Retinotopy'); % should be sym link in git folder

% Derive other file paths
dirPth.bs.dataPth = fullfile(dirPth.bsPth,'data',subjID);
dirPth.bs.anatPth = fullfile(dirPth.bsPth,'anat',subjID);

%% ----- Modelfitting ----

dirPth.model.saveFigPth   = fullfile(mprf_rootPath, 'data','Retinotopy', 'Quality_check', subjID, 'modelfit'); % should be sym link in git folder
dirPth.model.saveDataPth  = fullfile(mprf_rootPath, 'data','Retinotopy', 'Quality_check', subjID, 'modelfit');

end