function opt = getOpts(varargin)
% Function to get struct with default analysis pipeline options for MEG
% Retinotopy project. In case you want to change the default, use following
% example:
%
% Example 1:
%   opt = getOpts('verbose', false)
% Example 2:
%   opt = getOpts('foo', true)


% --- GENERAL ---
opt.skipMEGPreproc        = true;           % General
opt.skipMRIPreproc        = true;           % General
opt.verbose               = false;           % General
opt.doSaveData            = true;           % General
opt.saveFig               = true;           % General
opt.fullSizeGainMtx       = false;          % General: if true, execute analysis with fullsize meshes and gain matrix (FS size), if false, downsample to Brainstorm mesh size
opt.surfVisualize         = false;          % General: visualize surface meshes

% --- MEG Preproc ---
opt.doFiltering           = true;           % MEG preprocessing
opt.doDenoise             = true;           % MEG preprocessing
opt.removeStartOfRunEpoch = false;          % MEG preprocessing
opt.varThreshold          = [0.05 20];      % MEG preprocessing
opt.badChannelThreshold   = 0.2;            % MEG preprocessing
opt.badEpochThreshold     = 0.2;            % MEG preprocessing
opt.triggerChan           = 161:168;        % MEG Sensor info
opt.photoDiodeChan        = 192;            % MEG Sensor info
opt.dataChan              = 1:157;          % MEG Sensor info
opt.fs                    = 1000;           % MEG Sensor info: Sample rate (Hz)
opt.flickerFreq           = 10;             % MEG Experiment info: Stim freq (Hz)
opt.epochStartEnd         = [0.15 (0.15+1.1)]; % MEG Experiment info: Epoch length (s), first 150 ms are blank, one epoch length = 1.100 s,

% --- MRI pRF model ---
opt.eccThresh             = [0 10];         % MRI prf model: only use eccentricities that fall within the stimulus (20 deg in diameter)
opt.varExplThresh         = [0.1 inf];      % MRI prf model: Remove low variance explained vertices
opt.betaPrctileThresh     = [0 95];         % MRI prf model: Remove beta outliers by only taking up to the 95th percentile
opt.useSmoothedData       = true;           % MRI prf model: Use smoothed surface data or not
opt.useBensonMaps         = false;          % MRI prf model: Make prediction from Benson retinotopy atlas, instead of actual retinotopy data

% --- PERTUBATION of pRF models ---
opt.perturbOrigPRFs       = false;          % PERTUBATION of MRI prf model: say false for none, or choose from 'position', 'size', 'scramble'
if opt.perturbOrigPRFs
    opt.varyPosition          = -pi:(pi/4):pi;  % PERTUBATION of MRI prf model: Vary the prf center (polar angle rotation)
    opt.varySize              = unique(round(logspace(log10(0.2),log10(10),20),1)); % PERTUBATION of MRI prf model: Vary the prf size (scale sigma)
    opt.nScrambles            = 1000;
    opt.subfolder             = ['vary_' opt.perturbOrigPRFs]; % create subfolder in case of saving images / data
else
    opt.varyPosition          = [];
    opt.varySize              = [];
    opt.nScrambles            = [];
    opt.subfolder             = 'original'; 
end

% --- ROI (either draw manually on mrMesh and export to FS (FreeSurfer) and BS (BrainStorm) or Use wang et al atlas)
opt.roimrvToFS = 0; % set to 0 if rois drawn on FS surface are used

% --- Folders and filenames ---
opt.fNamePostFix          = sprintf('_benson%d_highres%d_smoothed%d', ...
    opt.useBensonMaps, opt.fullSizeGainMtx, opt.useSmoothedData);


%% Check for extra inputs in case changing the default options
if exist('varargin', 'var')
    
    % Get fieldnames
    fns = fieldnames(opt);
    for ii = 1:2:length(varargin)
        % paired parameter and value
        parname = varargin{ii};
        val     = varargin{ii+1};
        
        % check whether this parameter exists in the defaults
        idx = cellfind(fns, parname);
        
        % if so, replace it; if not add it to the end of opt
        if ~isempty(idx), opt.(fns{idx}) = val;
        else, opt.(parname) = val; end
        
        
    end
end

return