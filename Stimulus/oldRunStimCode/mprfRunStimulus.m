% To present the stimuli....
%
% To figure out: delay after pressing '5'. Screen keeps hanging for a while
% here and causes a delay in the stimulus sequence. As it seems, the code
% does not make up for this by speeding up the presentation of the first
% few frames. Is this because of the set up I am working on? Or of
% something else? -> moving to the winawer lab vistadisp solved it....
%
%

% To start things:
% addpath(genpath('/Users/bpklein/matlab/git/vistadisp'))
% tbUse('psychtoolbox-3');
% cd('/Volumes/server/Projects/MEG/Retinotopy');
% addpath(genpath(pwd));




type = 'meg_ff';




switch lower(type)
    
    
    case 'meg_ff'
        
        PTBInitStimTracker;
        global PTBTriggerLength
        PTBTriggerLength = 0.001;
        
        
        
        tr = 1.5;
        stimfile = 'MEG_FF_patterns_rand_perm_9.mat';
        
        cal = 'meg_lcd';
        d   = loadDisplayParams(cal);
        %         if any(size(d.gamma) == 1)
        %             d.gamma = reshape(d.gamma,length(d.gamma)/3,3);
        %         end
        
        
        params = retCreateDefaultGUIParams;         % Some default parameters
        % Check what these lines below do to the code...
        params.modality         = 'meg';
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
        
        
        params.loadMatrix = stimfile;
        ret(params)
        
        
    case 'fmri_carrier'
        
        
        tr = 1.5;
        stimfile = 'stimulus_fmri_carrier_dartboard_run1.mat';
        
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
        params.skipSyncTests    = 1;
        
        
        params.loadMatrix = stimfile;
        ret(params)
        
        
        
    case 'fmri_ret'
        
        
        tr = 1.5;
        stimfile = 'MRI_ret_stim_try.mat';
        
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
        params.skipSyncTests    = 1;
        
        
        params.loadMatrix = stimfile;
        ret(params)
        
        
        
        
        
end



