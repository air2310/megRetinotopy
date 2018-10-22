%% s_preprocessMRIRetinotopyData

% This script is the main analysis to preprocess the MRI dataset.

% This script relies on the following toolboxes:
% * Vistasoft


% 0. before this, you should have:
%    - ran FS recon-all on the T1
%    - Convert FS surfaces to Vista
%    - preprocess functional data


% 1. Init session
initVistaMegRet

% 2. ComputeMean
computeMeanMapCoronalMegRet

% 3. makeStimFiles
makeStimFiles

% 4. solve_PRFs
solvePRFsMegRet

% 5. Export Niftis to surface
exportNiftisMegRet


% 6. Run Bayesian template and Wang atlas with Noah Benson's docker, load
% these and then export ROIs from FS to Vista
...


