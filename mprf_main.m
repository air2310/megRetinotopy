function mprf_main(subjID, opt)
%
% Wrapper script containing MEG and MRI preprocessing and analyses subfunctions
% involved in the MEG Retinotopy project.
%
% %%% WORKFLOW %%%
% 0. Load paths and define parameters
%
% 1. MEG data preprocessing:
%   1.0: Define preprocessing options (if not defined as input variable)
%   1.1: Preprocess MEG data from raw (if requested)
%   1.2: Load MEG stim
%   1.3: Load MEG Gain matrix
%
% 2. MRI data preprocessing: (if requested)
%   2.1: Get pRF parameters in mrVista voxel space and smooth parameters
%   2.2: Project pRF parameters to freesurfer anatomy, get visual ROIs
%   2.3: Project pRF parameters from Freesurfer to Brainstorm surface
%   (2.4): if perturbing pRF parameters on surface, we do so here.
%
% 3. Forward model:
%   3.1: Predict response for MEG stimulus at cortical surface level
%           (could be Brainstorm or FreeSurfer surface)
%   3.2: Predict response for MEG stimulus at MEG sensor level
%           (multiply with gain matrix)
%   3.3: Computing phase referenced amplitude from preprocessed MEG data
%           and predicted MEG responses
%   3.4: Comparing predicted MEG time series and phase-referenced MEG 
%           steady-state responses
%
%
%
% %%% DEPENDENCIES %%%
% A. Preprocessing steps:
%   A1: FreeSurfer's auto-segmentation (v6?)
%   A2: preprocessing of raw MRI data into preprocessed nifti's (i.e. MRI
%       distortion correction w/ fsl top-up)
%   A3: Running pRF model in mrVista (voxel space) to get Gray Retinotopy modelfits
%   A4: Having Wang et al. (2015) atlas in subject's FreeSurfer directory
%
% B. External MATLAB Toolboxes:
% - Brainstorm (v??)
% - FieldTrip (v??)
% - VistaSoft (v??)
% - meg_utils (v??)
%
% Add all with the ToolboxToolbox:
%   tbUse('retmeg')
%
%
% Written by Akhil Edadan (UU) and Eline Kupers (NYU) - 2019
%

%% 0. Load paths

% Load paths with data files for this subject
dirPth = loadPaths(subjID);

% Go back to root
cd(mprf_rootPath)

% Set options if not defined (see getOpts function for more options)
if ~exist('opt', 'var') || isempty(opt,'var')
    opt = getOpts('saveFig', true,'verbose', true, 'fullSizeMesh', true, ...
        'perturbOrigPRFs', false, 'addOffsetParam', false, ...
        'refitGainParam', refitGain);
end

fprintf('(%s): Starting analysis of subject %s, using %s\n', mfilename, subjID, regexprep(opt.subfolder,'/',' '));

%% 1. MEG data preprocessing

% Start MEG data processing (if requested)
if ~opt.skipMEGPreproc
    if opt.verbose; fprintf('(%s): Preprocess MEG data..\n', mfilename); end
    
    % 1.1 Get preprocessed data from raw MEG data (.sqd) to preprocessed MEG data
    % (matfile, MEG sensors x epochs x time points x run nr)
    [data, conditions, opt] = preprocessMEGRetinotopyData(subjID, dirPth, opt);
    data = data.data;
    
    % 1.2 Get MEG stimulus (binarized and reduced to epochs x 10201 pixels)
    stim  = loadStim(subjID, dirPth, opt);
    stim  = stim.MEG;
    
    % 1.3 Get Gain matrix (produced via brainstorm)
    gainMtx  = loadGainMtx(subjID, dirPth, opt);
    
else % If you want to skip preprocessing, load data structs
    load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'epoched_data_hp_preproc_denoised.mat'), 'data');
    load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'megStimConditions.mat'), 'triggers');
    load(fullfile(dirPth.meg.processedDataPth, 'meg_stimulus.mat'), 'meg_stim');
    gainMtx    = loadGainMtx(subjID, dirPth, opt);
    
    % remove similar-named fields
    data       = data.data;
    stim       = meg_stim;
    conditions = triggers;
end

% Get MEG data struct and fill with loaded data
meg = struct();
meg.data            = data;
meg.stim            = stim;
meg.stim.conditions = conditions;
meg.gain            = gainMtx;

% Save some memory
clear gainMtx data stim conditions


%% 2. MRI data preprocessing
% This section will:
% 2.1 Smoothing pRF params in voxel space (and recompute beta's)
% 2.2. Export pRF params to FreeSurfer mesh + get Wang et al visual rois
% 2.3  Export pRF to Brainstorm mesh (if opt.fullSizeMesh = false)

% Get MRI data struct
mri     = struct();

% What pRF retinotopy maps are used on what size mesh?
if opt.mri.useBensonMaps % Use benson maps or not
    [mri.data, mri.prfSurfPath] = loadBensonRetinotopyMaps(subjID, dirPth, opt);
elseif opt.fullSizeMesh % if using full gain matrix, we need pRF params on FreeSurfer surface
    mri.prfSurfPath = dirPth.fmri.saveDataPth_prfFS;
else % if using downsampled gain matrix, we need pRF params on Brainstorm surface
    mri.prfSurfPath = dirPth.fmri.saveDataPth_prfBS;
end

% Start MRI data processing (if requested)
if ~opt.skipMRIPreproc
    if opt.verbose; fprintf('(%s): Preprocess MRI data..\n', mfilename); end
    
    % 2.1 Smoothing pRF params in voxel space (and then recompute beta 
    % values with smoothed parameters). We use mrVISTA for this.
    % mrVista gray nodes (i.e. voxels) --> mrVista gray nodes (i.e. voxels)
    mprf_pRF_sm(dirPth, opt);
    
    % Get summary figures for pRF parameters before/after smoothing
    if opt.verbose
        mprf_pRF_sm_fig(dirPth, opt);
        close all;
    end
    
    % 2.2 Get smoothed pRF params on FreeSurfer mid gray surface.
    % --> mrVista gray nodes (i.e. voxels) to FreeSurfer vertices
    mprf_pRF_sm_FS(dirPth,opt);
    
    % Get summary figures for pRF parameters on FreeSurfer surface
    if opt.verbose
        mprf_pRF_sm_FS_fig(dirPth,opt);
        close all;
    end
    
    % 2.3 Get smoothed pRF params and ROIs on Brainstorm pial surface.
    % --> Freesurfer vertices to downsampled Brainstorm vertices
    mprf_pRF_sm_FS_BS(subjID, dirPth,opt);
    
    % Get summary figures for pRF parameters on Brainstorm surface
    if opt.verbose
        mprf_pRF_sm_FS_BS_fig(dirPth,opt);
        close all;
    end
end


% 2.4 If perturbing original pRF parameters on the cortical surface:
if opt.vary.perturbOrigPRFs
    if opt.verbose; fprintf('(%s): Perturb local pRFs on cortex..\n', mfilename); end
    mprf_perturbOrigPRFs(mri.prfSurfPath, opt)
end


%% 3. Forward model

% 3.1 Predict response to MEG stimulus on mesh vertices (could be BS or FS)
%       INPUTS  (1) path to pRF parameters on surface (string)
%               (2) MEG stimulus (struct with x, y, images, etc)
%       OUTPUTS (1) predicted responses on surface (epochs x vertices)

predSurfResponse = mprf_MEGPredictionFromSurfaceWrapper(mri.prfSurfPath, meg.stim, dirPth, opt);

%% 3.2 Predicted response for MEG stimulus at MEG sensor level (weighting
%     predicted surface responses with gain matrix)
%       INPUTS  (1) predicted surface responses (on BS or FS surface)
%                                           (epochs x vertices)
%               (2) gain matrix             (sensors x vertices)
%
%       OUTPUTS (1) predicted MEG responses (epochs x sensors)

predMEGResponse = mprf_MEGPredictionSensorsWrapper(predSurfResponse, meg.gain, dirPth, opt);

%% 3.3 Computing phase referenced amplitude from preprocessed MEG data
% and predicted MEG responses from cortical surface
%   	INPUTS  (1) preprocessed MEG data   (time x epochs x run x sensors)
%               (2) predicted MEG responses (epochs x sensors)
%       OUTPUTS (1) phase-referenced MEG time series
%                                           (sensors x run group x epochs)
%               (2) bestBetas               (1 x run group x sensors)
%               (3) bestRefPhase            (1 x run group x sensors)
%               (4) bestOffsets             (1 x run group x sensors)


[phaseRefMEGResponse, bestBetas, bestRefPhase, bestOffsets] = mprf_MEGPhaseReferenceDataWrapper(meg.data, predMEGResponse, dirPth, opt);

% 3.4 Comparing predicted MEG time series and phase-referenced MEG steady-state responses
%       INPUTS  (1) Phase referenced MEG time series     (sensors x run groups x epochs)
%               (2) predicted MEG sensor responses from MRI prfs
%                                                        (epochs x sensors)
%       OUTPUTS (1) modelfit to mean phase-referenced MEG data, scaled by
%                   gain (and if requested, with offset) (epochs x sensors)
%               (2) average variance explained per MEG sensor (1 x sensors)

[predMEGResponseScaled,meanVarExpl] = mprf_CompareMEGDataToPRFPredictionWrapper(phaseRefMEGResponse, predMEGResponse, bestBetas, bestOffsets, dirPth, opt);


fprintf('(%s) finished without error!\n', mfilename)

end
