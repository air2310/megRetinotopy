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
function mprf_main(subjID)
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
% subjID = 'wlsubj111';

%%
% Load paths with data files for this subject
dirPth = loadPaths(subjID);

% Go back to root
cd(mprf_rootPath)

% Set options
opt = getOpts('saveFig',1,'verbose',1,'skipMRIPreproc',1,'skipMEGPreproc',1,'perturbOrigPRFs','position');
%opt = getOpts('saveFig',1,'verbose',1,'skipMRIPreproc',1,'skipMEGPreproc',0);
%opt = getOpts('saveFig',1,'verbose',1);

if strcmp(subjID, 'wlsubj039')
    opt = getOpts('useCoherentSpectrum', true','skipMRIPreproc',0,'verbose',0);
end

%% 1. MEG data preprocessing

if ~opt.skipMEGPreproc
    % 1.1 Get preprocessed data from raw MEG data (.sqd) to preprocessed MEG data
    % (matfile, MEG sensors x epochs x time points x run nr)
    [data, conditions, opt] = preprocessMEGRetinotopyData(subjID, dirPth, opt);
    data = data.data;
    
    % 1.2 Get MEG stimulus (binarized and reduced to epochs x 10201 pixels)
    stim  = loadStim(subjID, dirPth, opt);
    stim  = stim.MEG;
    
    % 1.3 Get Gain matrix (produced via brainstorm)
    gainMtx  = loadGainMtx(subjID, dirPth, opt);
    
else % If you want to skip preprocessing
    load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'epoched_data_hp_preproc_denoised.mat'), 'data');
    load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'megStimConditions.mat'), 'triggers');
    load(fullfile(dirPth.meg.processedDataPth, 'meg_stimulus.mat'), 'meg_stim');
    gainMtx    = loadGainMtx(subjID, dirPth, opt);
    
    % remove similar-named fields
    data       = data.data;
    stim       = meg_stim;
    conditions = triggers;
end

meg = struct();
meg.data            = data;
meg.stim            = stim;
meg.stim.conditions = conditions;
meg.gain            = gainMtx;

% Save some memory
clear gainMtx data stim conditions


%% 2 MRI data preprocessing

% Pre-requisits:
% (1) preprocessing of raw data into preprocessed nifti's
% (2) Running pRF model in mrVista (voxel space) to get Gray RM modelfits
% (3) Having Wang et al. (2015) atlas in subject's FreeSurfer directory 

mri     = struct();

if opt.mri.useBensonMaps % Use benson maps or not
    [mri.data, mri.prfSurfPath] = loadBensonRetinotopyMaps(subjID, dirPth, opt);
elseif opt.fullSizeMesh % if using full gain matrix, we need pRF params on FreeSurfer surface
    mri.prfSurfPath = dirPth.fmri.saveDataPth_prfFS;  
else % if using downsampled gain matrix, we need pRF params on Brainstorm surface
    mri.prfSurfPath = dirPth.fmri.saveDataPth_prfBS;
end

% This section will get Wang et al rois + smoothing pRF params (recomp
% beta) + exporting to FS + exporting to BS
% 
% input - pRF parameters in mrVista voxel space
%       - freesurfer directory for the subject (from step Anatomy)
%       - freesurfer anatomy - white and pial surface files
%       - BS pial surface files
% output - pRF parameters in BS

if ~opt.skipMRIPreproc

    % MRVISTA>>MRVISTA: Smoothing pRF params (and then recompute beta values with smoothed parameters)
    mprf_pRF_sm(dirPth, opt);

    % Get summary figures for pRF parameters before/after smoothing
    if opt.verbose
        mprf_pRF_sm_fig(dirPth, opt); 
        close all;
    end
    
    % MRVISTA>>FREESURFER: Get smoothed pRF params on FreeSurfer mid gray surface
    mprf_pRF_sm_FS(dirPth,opt);
    
    % Get summary figures for pRF parameters before/after smoothing
    if opt.verbose
        mprf_pRF_sm_FS_fig(dirPth,opt);
        close all;
    end
    
    % FREESURFER>>BRAINSTORM: Get smoothed pRF params and ROIs on Brainstorm pial surface
    mprf_pRF_sm_FS_BS(subjID, dirPth,opt);
    
    if opt.verbose
        mprf_pRF_sm_FS_BS_fig(dirPth,opt);
        close all;
    end
end


% If perturbing original pRF parameters on the cortical surface:
if opt.vary.perturbOrigPRFs  
    mprf_perturbOrigPRFs(mri.prfSurfPath, opt) 
end


%% 3. Forward model

% 3.1 Predict response to MEG stimulus on mesh vertices (could be BS or FS)
%       inputs (1) path to pRF parameters on surface (string)
%              (2) MEG stimulus (struct with x, y, images, etc)
%       output - predicted responses on surface (epochs x vertices)

predSurfResponse = mprf_MEGPredictionFromSurfaceWrapper(mri.prfSurfPath, meg.stim, dirPth, opt);

% 3.2 Predicted response for MEG stimulus at MEG sensor level (weighting
%     predicted surface responses with gain matrix)
%       inputs (1) predicted surface responses (on BS or FS surface)
%                   (epochs x vertices)
%              (2) gain matrix (sensors x vertices)
%       output - predicted MEG responses (epochs x sensors)

predMEGResponse = mprf_MEGPredictionSensorsWrapper(predSurfResponse, meg.gain, dirPth, opt);

% 3.3 Computing phase referenced amplitude from preprocessed MEG data 
% and predicted MEG responses from cortical surface
%   	inputs (1) preprocessed MEG data (time x epochs x run x sensors)
%              (2) predicted MEG responses (epochs x sensors)
%       output - Phase referenced MEG time series (sensors x epochs)

phaseRefMEGResponse = mprf_MEGPhaseReferenceDataWrapper(meg.data, predMEGResponse, dirPth, opt);

% 3.4 Comparing predicted MEG time series and phase-referenced MEG steady-state responses
%       inputs (1) Phase referenced MEG time series (sensors x time)
%              (2) predicted MEG sensor responses from MRI prfs
%       outputs(1) modelfit to mean phase-referenced MEG data,
%              (2) variance explained per MEG sensor 

[meanPredResponse,meanVarExpl] = mprf_CompareMEGDataToPRFPredictionWrapper(phaseRefMEGResponse, predMEGResponse, dirPth, opt);


fprintf('(%s) Done!', mfilename)

end
