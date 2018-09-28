function runme_MEG_PRF_Retinotopy_MRI_PRISMA(run)


tr = 1.5;
stimfile = sprintf('MRI_retinotopy_stimulus_run_%d.mat',run);

cal = 'CBI_Propixx';
d   = loadDisplayParams(cal);

% Hacked these lines from Eline's runme_MEG_OnOffLeftRight_ET_M2008

params = retCreateDefaultGUIParams;         % Some default parameters
% Check what these lines below do to the code...
params.modality         = 'fmri';
params.prescanDuration  = 0;
params.interleaves      = NaN;
params.tr               = tr;
params.calibration      = cal;
params.framePeriod      = tr;
params.startScan        = 0;
params.motionSteps      = 2;
params.tempFreq         = 6/tr;
params.repetitions      = 1;
params.experiment       = 'Experiment From File';
params.period           = 12*params.tr;
params.numCycles        = 6;
params.skipSyncTests    = 1;
params.triggerKey       = '5';
params.fixation         = 'dot with grid';

params.loadMatrix = stimfile;
ret(params)






end


















