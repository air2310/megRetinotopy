function mprf_addPaths()
% Add subfolder paths for MEG Retinotopy project

% Get location of file and project
currFilePath = which('mprf_addPaths');
projDir = fileparts(currFilePath);

% Go to that directory
cd(projDir)

% Add subfolders, except for data folder (since that contains a lot of
% other folders which are not necessary for this project and slow things
% down)
addpath(genpath('./analysis'))
addpath(genpath('./figurescripts'))
addpath(genpath('./Devel'))
addpath(genpath('./MEG_analysis'))
addpath(genpath('./mprfSession'))
addpath(genpath('./MRI_analysis'))
addpath(genpath('./preprocessing'))
addpath(genpath('./Server_analysis'))
addpath(genpath('./StimRefFwdMdl_analysis'))

return