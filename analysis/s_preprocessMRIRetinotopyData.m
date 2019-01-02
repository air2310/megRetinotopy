%% s_preprocessMRIRetinotopyData

% This script is the main analysis to preprocess the MRI dataset.

% This script relies on the following toolboxes:
% * Vistasoft

% 0. before this, you should have:
%    - ran FS recon-all on the T1
%    - create a t1_class.nii.gz file from the ribbon.mgz
%    - Convert FS T1 and surfaces to Vista
%    - preprocess functional data
%    - have the freesurfer subject directory defined as as an environmental
%    variable: setenv('SUBJECTS_DIR', '/Volumes/server/Freesurfer_subjects')

subject      = 'wlsubj040';
sessionDir   = sprintf('/Volumes/server/Projects/MEG/Retinotopy/Data/fMRI/%s/vistaSession/', subject);
sessionName  = 'raw';
bidsSession  = 'nyu3TAllegra'; %'nyu3t01' or 'nyu3TAllegra' for subj004 and subj040

%% 1. Init session
megRet_initVista(subject, sessionDir, bidsSession)

%% 2. ComputeMean
megRet_computeMeanMapCoronal(subject, sessionDir)

%% 3. makeStimFiles
megRet_makeStimFiles(subject, sessionDir, sessionName, bidsSession)

%% 4. solve_PRFs
megRet_solvePRFs(subject, sessionDir, sessionName)

%% 5. Export Niftis to surface
megRet_exportNiftis(subject, sessionDir)

%% 6. Run Bayesian template and Wang atlas with Noah Benson's docker, load
% these and then export ROIs from FS to Vista
megRet_getROIsFromTemplate(subject, fullfile(sessionDir))


