function runModelMain5(subjID)
%% ORIGINAL ANALYSIS - BENSON MAPS

%%
% Load paths with data files for this subject
dirPth = loadPaths(subjID);

% Go back to root
cd(mprf_rootPath)

% Set options
% opt = getOpts;

opt = getOpts('verbose',true, 'useBensonMaps', true', 'useSmoothedData', false);

if strcmp(subjID, 'wlsubj068')
    opt = getOpts('verbose',true, 'useBensonMaps', true', 'useSmoothedData', false, 'betaPrctileThresh', [0 94]); % default threshold of 95th percentile still gives very large/outlier beta values
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
    gainMtx             = loadGainMtx(subjID, dirPth, opt);
    
else % If you want to skip preprocessing
    data       = load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'epoched_data_hp_preproc_denoised.mat'));
    conditions = load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'megStimConditions.mat'));
    stim       = load(fullfile(dirPth.meg.processedDataPth, 'meg_stimulus.mat'));
    gainMtx    = loadGainMtx(subjID, dirPth, opt);
    
    % remove similar-named fields
    data       = data.data.data;
    stim       = stim.meg_stim;
    conditions = conditions.triggers;
end

meg = struct();
meg.data            = data;
meg.stim            = stim;
meg.stim.conditions = conditions;
meg.gain            = gainMtx;

% Save some memory
clear gainMtx data stim conditions


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
    mprf_pRF_sm(dirPth, opt); % pRF params from mrV >> smoothed pRF params on mrV (flag)

    if opt.verbose
        mprf_pRF_sm_fig(dirPth, opt); % Generates summary figures for the pRF parameters before after smoothing
    end
    
    mprf_pRF_sm_FS(dirPth,opt); % smoothed pRF params on mrV >> smoothed pRF params on FS
     
    if opt.verbose
        mprf_pRF_sm_FS_fig(dirPth,opt);
    end
    
    mprf_pRF_sm_FS_BS(dirPth,opt); % smoothed pRF params on FS >>  smoothed pRF params on BS
    
    if opt.verbose
        mprf_pRF_sm_FS_BS_fig(dirPth,opt);
    end
    % smoothed prf parameters + ROIs on BS (pial) surface saved

end

% If perturbing original pRF parameters on the cortical surface:
if opt.perturbOrigPRFs  
    mprf_perturbOrigPRFs(mri.prfSurfPath, opt) 
end


%% 3. Forward model

% 3.1 Predict response for MEG stimulus at Surface level (could be BS or FS)
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


%% Figures

% Figure 1. Time series (1a)
%           MEG head plot (1b)

%Figure 2. Position range line plot
%          headplots for every position range
% makeFigure2(subjID, dirPth, opt);

%Figure 3. Size range line plot
%          headplots for every size range
% makeFigure3(subjID, dirPth, opt);



end