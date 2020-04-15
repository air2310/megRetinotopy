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
opt.skipMEGPreproc        = true;               % General
opt.skipMRIPreproc        = true;               % General
opt.verbose               = true;               % General
opt.saveData              = true;               % General
opt.saveFig               = true;               % General
opt.fullSizeMesh          = false;              % General: if true, execute analysis with fullsize meshes and gain matrix (FS size), if false, downsample to Brainstorm mesh size
opt.surfVisualize         = false;              % General: visualize surface meshes
opt.subfolder             = 'original';         % General: Create subfolder to save figures
opt.makeAverageFig        = false;              % General: Plot average across subjects
opt.headmodel             = 'OS';               % General: choose from 'OS' (default) or 'BEM' head model (for now, BEM is only implemented for wlsubj081 and wlsubj111)
opt.addOffsetParam        = true;               % General: if true, use both gain and offset parameters in fit, if false, regress with only 1 free param (gain)
opt.refitGainParam        = true;               % General: if true, recompute the gain factor for the prediction by refitting the phase-ref data.
opt.seed                  = 1;                  % General: random number generator seed for reproducing results, can be an integer > 0  or 'shuffle' to get reset rng

% --- MEG Preproc ---
opt.meg.doFiltering           = false;          % MEG preprocessing
opt.meg.doDenoise             = true;           % MEG preprocessing
opt.meg.removeStartOfRunEpoch = false;          % MEG preprocessing
opt.meg.varThreshold          = [0.05 20];      % MEG preprocessing
opt.meg.badChannelThreshold   = 0.2;            % MEG preprocessing
opt.meg.badEpochThreshold     = 0.2;            % MEG preprocessing
opt.meg.triggerChan           = 161:168;        % MEG Sensor info
opt.meg.photoDiodeChan        = 192;            % MEG Sensor info
opt.meg.dataChan              = 1:157;          % MEG Sensor info
opt.meg.fs                    = 1000;           % MEG Sensor info: Sample rate (Hz)
opt.meg.flickerFreq           = 10;             % MEG Experiment info: Stim freq (Hz)
opt.meg.epochStartEnd         = [0.15 (0.15+1.1)]; % MEG Experiment info: Epoch length (s), first 150 ms are blank, one epoch length = 1.100 s,
opt.meg.useCoherentSpectrum   = true;           % MEG Data for prf model: Use the coherent spectrum (average before FFT), instead of the incoherent spectrum (average after FFT)
opt.meg.doSplitHalfReliability= false;          % MEG amplitude split half correlation (1000x computed).

% --- MRI pRF model ---
opt.mri.eccThresh             = [0 10];         % MRI prf model: only use eccentricities that fall within the stimulus (20 deg in diameter)
opt.mri.varExplThresh         = [0.1 inf];      % MRI prf model: Remove low variance explained vertices
opt.mri.betaPrctileThresh     = [0 100];        % MRI prf model: Remove beta outliers by only taking data between 0 and xxth percentile
opt.mri.useSmoothedData       = true;           % MRI prf model: Use smoothed surface data or not
opt.mri.useBensonMaps         = false;          % MRI prf model: Make prediction from Benson retinotopy atlas, instead of actual retinotopy data
opt.mri.predSurfVarThresh     = [0 500];        % MRI prf model: if the variance of a predicted vertex response response is x1 times smaller or x2 times bigger than median, remove this response 

% --- PERTUBATION of pRF models ---
opt.vary.perturbOrigPRFs       = false;         % PERTUBATION of MRI prf model: say false for none, or choose from 'position', 'size', 'scramble'
opt.vary.position              = [];
opt.vary.size                  = [];
opt.vary.nScrambles            = [];
    
% --- ROI (either draw manually on mrMesh and export to FS (FreeSurfer) and BS (BrainStorm) or Use wang et al atlas)
opt.roi.roimrvToFS            = false; % set to false, use wang atlas, otherwise rois manually drawn in mrVista
opt.roi.onlyV123WangAtlas     = false; % set to true in case you want to only use V1-V3 ROIs of the Wang atlas (for comparing against benson atlas results)

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
        
        elseif isempty(idx)
        
            % check whether this parameter exists in the substructs
            for ii = 1:length(fns)
                if isstruct(opt.(fns{ii}))
                    subfns = fieldnames(opt.(fns{ii}));
                    idx = cellfind(subfns, parname);
                    if ~isempty(idx)
                        opt.(fns{ii}).(subfns{idx}) = val;
                    end
                end
            end
        
        elseif isempty(idx)
            opt.(parname) = val;
        end  
                    
    end
end

% check in case option to perturb original prfs was turned to true again
if opt.vary.perturbOrigPRFs
    opt.vary.position          = -pi:(pi/4):pi;  % PERTUBATION of MRI prf model: Vary the prf center (polar angle rotation)
    opt.vary.size              = unique(round(logspace(log10(0.2),log10(10),20),1)); % PERTUBATION of MRI prf model: Vary the prf size (scale sigma)
    opt.vary.nScrambles        = 1000;
    opt.subfolder              = ['vary_' opt.vary.perturbOrigPRFs]; % create subfolder in case of saving images / data
end

% --- Folders and filenames ---
opt.fNamePostFix          = sprintf('_benson%d_highres%d_smoothed%d_%s', ...
                            opt.mri.useBensonMaps, opt.fullSizeMesh, opt.mri.useSmoothedData, opt.headmodel);


% check in case additional original analyses are performed, we want to save those results separately
if opt.roi.onlyV123WangAtlas
    opt.subfolder = 'original/onlyV123WangAtlas'; % add results to extra folder to not override original
elseif opt.mri.useBensonMaps
    opt.subfolder = 'original/BensonMaps'; % add results to extra folder to not override 
end

if opt.meg.useCoherentSpectrum
    opt.subfolder = [opt.subfolder '/coherent'];
else
    opt.subfolder = [opt.subfolder '/incoherent'];
end

if opt.addOffsetParam
    opt.subfolder = [opt.subfolder 'WithOffset'];
else
    opt.subfolder = [opt.subfolder 'NoOffset'];
end

if opt.refitGainParam
    opt.subfolder = [opt.subfolder 'RecompFinalBeta'];
end

if strcmp(opt.headmodel, 'BEM')
    opt.subfolder = [opt.subfolder '/BEM'];
end

if opt.fullSizeMesh
    opt.subfolder = [opt.subfolder '/FSMesh'];
end

return