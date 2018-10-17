%% s_preprocessMRIRetinotopyData

% This script is the main analysis to preprocess the MRI dataset.

% This script relies on the following toolboxes:
% * Vistasoft


% 0. before this, you should have:
%    - ran FS recon-all on the T1
%    - Convert FS surfaces to Vista
%    - preprocess functional data


% 1. Init session

% 2. ComputeMean

% 3. makeStimFiles

% 4. solve_PRFs

% 5. Export Niftis to surface

% 6. Run Bayesian template and Wang atlas

% 7. Load ROIs and export to Vista


