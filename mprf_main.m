%% mprf_main
%
% Wrapper script containing MEG and MRI preprocessing and analyses subfunctions
% involved in the MEG Retinotopy project.
%
% WORKFLOW:
% 0. Load paths and define parameters
%
% 1. MEG data preprocessing:
%   1.0 Define preprocessing options
%   1.1 Preprocess MEG data from raw
%   1.2 Load MEG stim
%   1.3 Load MEG Gain matrix
%
% 2. MRI data preprocessing
%   2.1: xx
%
% 3. Forward model
%   3.1: Predict response for MEG stimulus at Surface level (could be BS or FS)
%   3.2: Predict response for MEG stimulus at MEG sensor level (multiply
%           with gain matrix)
%   3.3: Computing phase referenced amplitude from preprocessed MEG data 
%           and predicted MEG responses from cortical surface
%
% DEPENDENCIES:
% 1. Preprocessing steps:
% - FreeSurfer's auto-segmentation (v6???)
% - MRI distortion correction w/ top-up??
%
% 2. Toolboxes:
% - Brainstorm (v??)
% - FieldTrip (v??)
% - VistaSoft (v??)
% - meg_utils (v??)
%
% Add all with the ToolboxToolbox:
%   tbUse('retmeg')
%
% 
% By Akhil Edadan (UU) and Eline Kupers (NYU) - 2019
%
%% 0. Load paths

% Define subject ID
subjID = 'wlsubj058';

%%
% Load paths with data files for this subject
dirPth = loadPaths(subjID);

% Go back to root
cd(mprf_rootPath)

% Set options
opt.skipMEGPreproc        = true;  % General
opt.skipMRIPreproc        = true;  % General
opt.verbose               = true;  % General
opt.doSaveData            = true;  % General
opt.saveFig               = true; % General
opt.fullSizeGainMtx       = true; % General

opt.doFiltering           = true;       % MEG preprocessing
opt.doDenoise             = true;       % MEG preprocessing
opt.removeStartOfRunEpoch = false;      % MEG preprocessing
opt.varThreshold          = [0.05 20];  % MEG preprocessing
opt.badChannelThreshold   = 0.2;        % MEG preprocessing
opt.badEpochThreshold     = 0.2;        % MEG preprocessing
opt.triggerChan           = 161:168;    % MEG Sensor info
opt.photoDiodeChan        = 192;        % MEG Sensor info          
opt.dataChan              = 1:157;      % MEG Sensor info
opt.fs                    = 1000;       % MEG Sensor info: Sample rate (Hz) 
opt.flickerFreq           = 10;         % MEG Experiment info: Stim freq (Hz)   
opt.epochStartEnd         = [0.15 (0.15+1.1)]; % MEG Experiment info: Epoch length (s), first 150 ms are blank, one epoch length = 1.100 s,

opt.eccThresh             = [0 10];    % MRI prf model: only use eccentricities that fall within the stimulus (20 deg in diameter)
opt.varExplThresh         = [0.1 inf]; % MRI prf model: Remove low variance explained vertices
opt.betaPrctileThresh     = [0 95];    % MRI prf model: Remove beta outliers by only taking up to the 95th percentile
opt.useSmoothedData       = true;      % MRI prf model: Use smoothed surface data or not

opt.useBensonMaps         = false;     % MRI prf model: Make prediction from Benson retinotopy atlas, instead of actual retinotopy data

%% 1. MEG data preprocessing

if ~opt.skipMEGPreproc
    % 1.1 Get preprocessed data from raw MEG data (.sqd) to preprocessed MEG data
    % (matfile, MEG sensors x epochs x time points x run nr)
    [data, conditions, opt] = preprocessMEGRetinotopyData(subjID, dirPth, opt);
    
    % 1.2 Get MEG stimulus (binarized and reduced to epochs x 10201 pixels)
    stim  = loadStim(subjID, dirPth, opt);
    
    % 1.3 Get Gain matrix (produced via brainstorm)
    gainMtx             = loadGainMtx(subjID, dirPth, opt);
    
else % If you want to skip preprocessing
    data       = load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'epoched_data_hp_preproc_denoised.mat'));
    conditions = load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'megStimConditions.mat'));
    stim       = load(fullfile(dirPth.meg.processedDataPth, 'meg_stimulus.mat'));
    gainMtx    = loadGainMtx(subjID, dirPth, opt);
end

meg      = struct();
meg.data = data.data.data;
meg.stim = stim;
meg.gain = gainMtx;


%% 2 MRI data preprocessing

mri     = struct();

if opt.useBensonMaps 
    [mri.data, mri.prfSurfPath] = loadBensonRetinotopyMaps(subjID, dirPth, opt);
elseif opt.fullSizeGainMtx
    mri.prfSurfPath = dirPth.fmri.saveDataPth_prfFS;  
else
    mri.prfSurfPath = dirPth.fmri.saveDataPth_prfBS;
end

% Preprocessing
% input - raw fMRI data (dicoms)
% output - preprocessed fMRI (.nii) - voxels x #tr

%************<code>***************

% Running pRF model in mrVista (voxel space)
% input - preprocessed fMRI (.nii)
%       - fMRI stimuli
% output - pRF parameters - voxels x #parameters

%************<code>***************

% Wang et al rois + smoothing pRF params (recomp beta) + exporting to FS + exporting to BS
% input - pRF parameters in mrVista voxel space
%       - freesurfer directory for the subject (from step Anatomy)
%         - anatomy
%         - pial surface files
%       - BS surface files
% output - pRF parameters in BS

%************<code>***************

% prf parameters + ROIs on mrVista volume space

if ~opt.skipMRIPreproc

    %mprf_ROI % ROIs on mrVista space
    mprf_pRF_sm(dirPth,opt) % pRF params from mrV >> smoothed pRF params on mrV (flag)

    mprf_pRF_sm_fig(subjID, dirPth, opt); % Generates summary figures for the pRF parameters before after smoothing

    mprf_pRF_sm_FS(dirPth) % smoothed pRF params on mrV >> smoothed pRF params on FS
    mprf_pRF_sm_FS_fig(dirPth);

    mprf_pRF_sm_FS_BS(dirPth) % smoothed pRF params on FS >>  smoothed pRF params on BS
    mprf_pRF_sm_FS_BS_fig(dirPth);

    % smoothed prf parameters + ROIs on BS (pial) surface saved

end
%% 3. Forward model

% 3.1 Predict response for MEG stimulus at Surface level (could be BS or FS)
%       inputs (1) path to pRF parameters on surface (string)
%              (2) MEG stimulus (struct with x, y, images, etc)
%       output - predicted responses on surface (epochs x vertices)

predSurfResponse = mprf_MEGPredictionFromSurface(mri.prfSurfPath, meg.stim, dirPth, opt);

% 3.2 Predicted response for MEG stimulus at MEG sensor level (weighting
%     predicted surface responses with gain matrix)
%       inputs (1) predicted surface responses (on BS or FS surface)
%                   (epochs x vertices)
%              (2) gain matrix (sensors x vertices)
%       output - predicted MEG responses (epochs x sensors)

predMEGResponse = mprf_MEGPredictionSensors(predSurfResponse, meg.gain, dirPth, opt);

% 3.3 Computing phase referenced amplitude from preprocessed MEG data 
% and predicted MEG responses from cortical surface
%   	inputs (1) preprocessed MEG data (time x epochs x run x sensors)
%              (2) predicted MEG responses (epochs x sensors)
%       output - Phase referenced MEG time series (sensors x epochs)

phaseRefMEGResponse = mprf_MEGPhaseReferenceData(meg.data, predMEGResponse, dirPth, opt);

% 3.4 Comparing predicted MEG time series and phase-referenced MEG steady-state responses
%       inputs (1) Phase referenced MEG time series (sensors x time)
%              (2) predicted MEG sensor responses from MRI prfs
%       outputs(1) modelfit to mean phase-referenced MEG data,
%              (2) variance explained per MEG sensor 

[meanPredResponse,meanVarExpl] = mprf_CompareMEGDataToPredictionFromMRIPRFs(phaseRefMEGResponse, predMEGResponse, dirPth, opt);


%% Figures

% Figure 1. Time series (1a)
%           MEG head plot (1b)

%Figure 2. Position range line plot
%          headplots for every position range

%Figure 2. Size range line plot
%          headplots for every size range


