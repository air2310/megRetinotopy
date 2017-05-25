function runme_2017_03_14_fMRI_carrier(run, carrier)

% Lenght of experiment 156 seconds with TR = 1.5 equals  104 TR
tr = 1.5;
stimfile = sprintf('stimulus_fmri_carrier_%s_run%d.mat',carrier,run);


cal = 'CBI_NYU_projector';
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
params.skipSyncTests    = 0;
params.triggerKey       = '`';
params.fixation         = 'dot with grid';


params.loadMatrix = stimfile;
ret(params)


end


