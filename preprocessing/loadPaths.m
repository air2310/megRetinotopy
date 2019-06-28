function dirPth = loadPaths(subjID)

dirPth  = struct();
dirPth.rootPth = pwd;

dirPth.subjID = subjID;


%% ------ MEG ------
dirPth.meg.dataPth   = './Data/MEG/';
dirPth.meg.saveFigPth    = fullfile('./Quality_check/',subjID, 'meg');

% Derive other file paths
dirPth.meg.rawSqdPth = fullfile(dirPth.meg.dataPth, subjID, 'raw');
dirPth.meg.paramFilePth = fullfile(dirPth.meg.dataPth, subjID,  'paramFiles');
dirPth.meg.stimFilePth = fullfile(dirPth.meg.dataPth, subjID, 'stimFiles');


%% ------ Anatomy ------ 
dirPth.mri.dataPth = './Data/Anatomy/';

% Derive other file paths
dirPth.mri.anatPth = fullfile(dirPth.mri.dataPth, subjID);
dirPth.mri.segPth = fullfile(dirPth.mri.dataPth, subjID);

%% ------ fMRI ------ 
dirPth.fmri.dataPth = './Data/fMRI/';
dirPth.fmri.saveDataPth    = fullfile('./Quality_check/',subjID, 'fmri');

% Derive other file paths
dirPth.fmri.mrvPth = fullfile(dirPth.fmri.dataPth, subjID, 'vistaSession');
dirPth.fmri.paramFilePth = fullfile(dirPth.fmri.dataPth, subjID,  'paramFiles');
dirPth.fmri.stimFilePth = fullfile(dirPth.fmri.dataPth, subjID, 'stimFiles');

dirPth.fmri.saveDataPth_prfMrv = fullfile(dirPth.fmri.saveDataPth, 'prfMrv');
dirPth.fmri.saveDataPth_prfFS = fullfile(dirPth.fmri.saveDataPth, 'prfFS');
dirPth.fmri.saveDataPth_prfBS = fullfile(dirPth.fmri.saveDataPth, 'prfBS');

dirPth.fmri.saveDataPth_roiFS = fullfile(dirPth.fmri.saveDataPth, 'roiFS');
dirPth.fmri.saveDataPth_roiBS = fullfile(dirPth.fmri.saveDataPth, 'roiBS');

%% ------ freesurfer ------ 
dirPth.freesurfer.dataPth = './Data/Freesurfer_directory/Freesurfer_subjects/';

% Derive other file paths
dirPth.freesurfer.surfPth = fullfile(dirPth.freesurfer.dataPth, subjID, 'surf');


%% ------ brainstorm ------ 
dirPth.brainstorm.dbPth = './Data/Brainstorm_db/';

% Derive other file paths
dirPth.brainstorm.dataPth = fullfile(dirPth.brainstorm.dbPth,'data',subjID);
dirPth.brainstorm.anatPth = fullfile(dirPth.brainstorm.dbPth,'anat',subjID);

