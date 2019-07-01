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
%   3.1: xx
%
%
% This script assumes the following preprocessing steps:
% - FreeSurfer's auto-segmentation (v6???)
% - MRI distortion correction w/ top-up??
%
% This script requires the following toolboxes:
% - Brainstorm (v??)
% - FieldTrip (v??)
% - VistaSoft (v??)
% - meg_utils (v??)
%
% Add all with the ToolboxToolbox:
%   tbUse('retmeg') 
%
% By Akhil Edadan (UU) and Eline Kupers (NYU) - 2019
%
%
%% 0. Load paths

% Define subject ID
subjID = 'wlsubj058';

% Load paths with data files for this subject
dirPth = loadPaths(subjID);
   

%% MEG 

% 1.0 Preprocessing options:
opt.verbose               = true;
opt.doFiltering           = true;
opt.doDenoise             = true;
opt.doSaveData            = true;
opt.saveFig               = false; 
opt.removeStartOfRunEpoch = false;

% 1.1 Get preprocessed data from raw MEG data (.sqd) to preprocessed MEG data
% (matfile, MEG sensors x epochs x time points x run nr) 
data        = preprocessMEGRetinotopyData(subjID, dirPth, opt);

% 1.2 Get MEG stimulus (binarized and reduced to epochs x 10201 pixels)
stim        = loadStim(subjID, dirPth, opt);

% 1.3 Get Gain matrix (produced via brainstorm)
opt.fullSizeGainMtx = false;
gainMtx             = loadGainMtx(subjID, dirPth, opt);

meg = struct();
meg.data = data;
meg.stim = stim;
meg.gain = gainMtx;

   

  
   
   
     
%% fMRI

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
     
     
     %mprf_ROI % ROIs on mrVista space
     plot_stim = 0; % flag to plot the stimulus movie
     mprf_pRF_sm(dirPth,plot_stim) % pRF params from mrV >> smoothed pRF params on mrV (flag)
     mprf_pRF_sm_fig(dirPth); % Generates summary figures for the pRF parameters before after smoothing
     
     mprf_pRF_sm_FS(dirPth) % smoothed pRF params on mrV >> smoothed pRF params on FS
     mprf_pRF_sm_FS_fig(dirPth);
     
     mprf_pRF_sm_FS_BS(dirPth) % smoothed pRF params on FS >>  smoothed pRF params on BS
     mprf_pRF_sm_FS_BS_fig(dirPth);
     
     % smoothed prf parameters + ROIs on BS (pial) surface saved 
     
%% Stimulus referred forward model

   % Predicted MEG response on BS  
     % input - pRF parameters in BS
     %       - MEG stimuli 
     % output - predicted MEG time series on BS
     
     %************<code>***************
     
     [prf,model] = mprf_MEG_pred_BS(subjID);
     
     % predicted response to MEG stimulus at every BS surface saved
     
   % Weighting response with gain matrix
     % input - predicted MEG time series on BS
     %       - gain matrix
     % output - predicted MEG time series on MEG sensor space
     
     %************<code>***************
     
     meg_resp = mprf_MEG_pred_Msen(subjID); 
     
     
   % Computing phase referenced amplitude from preprocessed MEG time series and
   % predicted MEG time series on MEG sensor space
     % input - preprocessed MEG (.m) - sensors x epochs x time x run 
     %       - predicted MEG time series on MEG sensor space
     % output - Phase referenced MEG time series (sensors x time)
   
     %************<code>***************
     
     pred.prf = prf;
     pred.model = model;
     pred.syn = 0;
     pred.meg_resp = meg_resp;

     cur_time = datestr(now);
     cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
     pred.cur_time = cur_time;
     
     mprfSession_run_original_model_server(pred,meg_data_file_path,plot_figure);

     
     mprf_PhRefAmp
     
   % Comparing predicted MEG time series and phase referenced MEG time series  
     % input - predicted MEG time series on MEG sensor space
     %       - Phase referenced MEG time series (sensors x time)   
     % output - variance explained per MEG sensor
   
     %************<code>***************
     
     
%% Figures
  
  % Figure 1. Time series (1a)
  %           MEG head plot (1b)
  
  %Figure 2. Position range line plot
  %          headplots for every position range          
  
  %Figure 2. Size range line plot
  %          headplots for every size range           
   
   
  