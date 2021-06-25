function dirPth = loadPaths(subjID)
%
% Function to create a struct with all different directories for the
% different types of files containing data

%% ----- General -----

dirPth  = struct();
dirPth.rootPth      = mprf_rootPath;
dirPth.subjID       = subjID;
retinoFolder        = fullfile(mprf_rootPath,'data','Retinotopy-1'); % should be sym link in git folder
dirPth.fsPth        = fullfile(mprf_rootPath,'data','Freesurfer_subjects-1'); % should be sym link in git folder
dirPth.bsPth        = fullfile(mprf_rootPath,'data','brainstorm_db-1','MEG_Retinotopy'); % should be sym link in git folder

%% Derive other paths
dirPth.sessionPth       = fullfile(retinoFolder, 'Subject_sessions', subjID);

%% ------ MEG ------
dirPth.meg.dataPth      = fullfile(retinoFolder, 'Data', 'MEG'); % should be sym link in git folder
dirPth.meg.saveFigPth   = fullfile(retinoFolder, 'Modelfits', subjID, 'meg'); % should be sym link in git folder

% Derive other file paths
dirPth.meg.rawSqdPth    = fullfile(dirPth.meg.dataPth, subjID, 'raw');
dirPth.meg.paramFilePth = fullfile(dirPth.meg.dataPth, subjID,  'paramFiles');
dirPth.meg.stimFilePth  = fullfile(dirPth.meg.dataPth, subjID, 'stimFiles');
dirPth.meg.processedDataPth = fullfile(dirPth.meg.dataPth, subjID, 'processed');

dirPth.meg.stimFile     = fullfile(dirPth.meg.stimFilePth, 'MEG_retinotopy_stimulus_run_1.mat');
dirPth.meg.stimGridFile = fullfile(dirPth.meg.stimFilePth, 'MEG_grid.mat');

dirPth.meg.eyePth       = fullfile(dirPth.meg.dataPth, subjID, 'eye');
%% ------ FreeSurfer ------ 
if strcmp(subjID, 'wlsubj039') || strcmp(subjID, 'wlsubj081')
    dirPth.fs.segPth    = fullfile(dirPth.fsPth, [subjID '_wVitaminE']);
else
    dirPth.fs.segPth    = fullfile(dirPth.fsPth, subjID);
end
dirPth.fs.surfPth       = fullfile(dirPth.fs.segPth, 'surf');


%% ------ fMRI ------ 
dirPth.fmri.dataPth     = fullfile(retinoFolder, 'Data', 'fMRI'); % should be sym link in git folder
dirPth.fmri.saveDataPth = fullfile(retinoFolder, 'Modelfits', subjID, 'fmri'); % should be sym link in git folder

dirPth.fmri.mrvPth       = fullfile(dirPth.fmri.dataPth, subjID, 'vistaSession');
dirPth.fmri.paramFilePth = fullfile(dirPth.fmri.mrvPth, 'Stimuli', 'paramFiles');
dirPth.fmri.stimFilePth  = fullfile(dirPth.fmri.mrvPth, 'Stimuli', 'stimFiles');
dirPth.fmri.eyePth       = fullfile(dirPth.fmri.mrvPth, 'eye');

dirPth.fmri.saveDataPth_prfMrv = fullfile(dirPth.fmri.saveDataPth, 'prfMrv');
dirPth.fmri.saveDataPth_prfFS = fullfile(dirPth.fmri.saveDataPth, 'prfFS');
dirPth.fmri.saveDataPth_prfBS = fullfile(dirPth.fmri.saveDataPth, 'prfBS');

dirPth.fmri.saveDataPth_roiFS = fullfile(dirPth.fmri.saveDataPth, 'roiFS');
dirPth.fmri.saveDataPth_roiBS = fullfile(dirPth.fmri.saveDataPth, 'roiBS');

d = dir(fullfile(dirPth.fmri.mrvPth, 'Gray', 'Averages','*fFit*'));
dirPth.fmri.vistaGrayFitFile  = fullfile(d.folder, d.name);

%% ------ Brainstorm ------ 

dirPth.bs.dataPth = fullfile(dirPth.bsPth,'data',subjID);
dirPth.bs.anatPth = fullfile(dirPth.bsPth,'anat',subjID);

%% ----- Modelfitting ----

dirPth.model.saveFigPth   = fullfile(retinoFolder, 'Modelfits', subjID, 'modelfit'); % should be sym link in git folder
dirPth.model.saveDataPth  = fullfile(retinoFolder, 'Modelfits', subjID, 'modelfit');

%% ----- Final figures ----

dirPth.finalFig.savePth = fullfile(retinoFolder, 'Modelfits', subjID, 'finalfig'); % should be sym link in git folder
dirPth.finalFig.savePthAverage = fullfile(retinoFolder, 'Modelfits', 'average','finalfig'); % should be sym link in git folder


end