%% s_preprocessMRIRetinotopyData

% This script is the main analysis to preprocess the MRI dataset.

% This script relies on the following toolboxes:
% * Vistasoft


% 0. before this, you should have:
%    - ran FS recon-all on the T1
%    - Convert FS surfaces to Vista
%    - preprocess functional data

subject      = 'wlsubj030';
sessionDir   = '/Volumes/server/Projects/MEG/Retinotopy/Data/fMRI/wlsubj030';
sessionName  = 'MRI_Data';
bidsSession  = 'nyu3T01';

% 1. Init session
<<<<<<< HEAD
megRet_initVista(subject)

% 2. ComputeMean
megRet_computeMeanMapCoronal(subject, sessionDir)

% 3. makeStimFiles
megRet_makeStimFiles(subject, sessionDir, sessionName, bidsSession)

% 4. solve_PRFs
megRet_solvePRFs(subject, sessionDir, sessionName, bidsSession)

% 5. Export Niftis to surface
exportNiftisMegRet

% 6. Run Bayesian template and Wang atlas with Noah Benson's docker, load
% these and then export ROIs from FS to Vista



