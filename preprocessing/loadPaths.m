function dirPth = loadPaths(subjID)
%
% Function to create a struct with all different directories for the
% different types of files containing data

%% ----- General -----
dirPth  = struct();
dirPth.rootPth = mprf_rootPath;
dirPth.subjID  = subjID;


%% ------ MEG ------
dirPth.meg.dataPth      = fullfile(mprf_rootPath, 'data','Retinotopy', 'Data', 'MEG');
dirPth.meg.saveFigPth   = fullfile(mprf_rootPath, 'data','Retinotopy', 'Quality_check', subjID, 'meg');

% Derive other file paths
dirPth.meg.rawSqdPth    = fullfile(dirPth.meg.dataPth, subjID, 'raw');
dirPth.meg.paramFilePth = fullfile(dirPth.meg.dataPth, subjID,  'paramFiles');
dirPth.meg.stimFilePth  = fullfile(dirPth.meg.dataPth, subjID, 'stimFiles');


%% ------ FreeSurfer ------ 
dirPth.fsPth            = fullfile(mprf_rootPath,'data','Freesurfer_subjects');

% Derive other file paths
dirPth.fs.segPth        = fullfile(dirPth.fsPth, subjID);
dirPth.fs.surfPth       = fullfile(dirPth.fs.segPth, 'surf');


%% ------ fMRI ------ 
dirPth.fmri.dataPth     = fullfile(mprf_rootPath, 'data','Retinotopy', 'Data', 'fMRI');
dirPth.fmri.saveDataPth = fullfile(mprf_rootPath, 'data','Retinotopy', 'Quality_check', subjID, 'fmri');

% Derive other file paths
dirPth.fmri.mrvPth       = fullfile(dirPth.fmri.dataPth, subjID, 'vistaSession');
dirPth.fmri.paramFilePth = fullfile(dirPth.fmri.mrvPth, 'Stimuli', 'paramFiles');
dirPth.fmri.stimFilePth  = fullfile(dirPth.fmri.mrvPth, 'Stimuli', 'stimFiles');

dirPth.fmri.saveDataPth_prfMrv = fullfile(dirPth.fmri.saveDataPth, 'prfMrv');
dirPth.fmri.saveDataPth_prfFS = fullfile(dirPth.fmri.saveDataPth, 'prfFS');
dirPth.fmri.saveDataPth_prfBS = fullfile(dirPth.fmri.saveDataPth, 'prfBS');

dirPth.fmri.saveDataPth_roiFS = fullfile(dirPth.fmri.saveDataPth, 'roiFS');
dirPth.fmri.saveDataPth_roiBS = fullfile(dirPth.fmri.saveDataPth, 'roiBS');

%% ------ Brainstorm ------ 
dirPth.bsPth = fullfile(mprf_rootPath, 'data', 'brainstorm_db', 'MEG_Retinotopy');

% Derive other file paths
dirPth.bs.dataPth = fullfile(dirPth.bsPth,'data',subjID);
dirPth.bs.anatPth = fullfile(dirPth.bsPth,'anat',subjID);

