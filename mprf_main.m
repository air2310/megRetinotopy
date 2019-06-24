% Script containing basic steps involved in the model

%% Anatomy
   % Segmentation (freesurfer ???)


%% MEG 

   % Preprocessing
     % input - raw MEG data (.sqd)
     % output - preprocessed MEG (.m) - sensors x epochs x time x run 

     %************<code>***************
  
   % MEG stimulus (binarized and reduced to epochs x 10201 pixels)
   
   % Gain matrix (from brainstorm)
     
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
     
     subjid='wlsubj004';
     
     %mprf_ROI % ROIs on mrVista space
     
     dir_pth.Anat_dir = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Data/Anatomy/%s',subjid);
     dir_pth.mrSession_dir = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Data/fMRI/%s/vistaSession',subjid);
     dir_pth.prf_dir = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Quality_check/%s/prf_data',subjid); % folder to save prf params (contains folders- nifti(.nii files), data_dir(.mat files))
     
     mprf_pRF_sm(subjid,dir_pth,1) % pRF params from mrV >> smoothed pRF params on mrV

     mprf_pRF_sm_FS(subjid) % smoothed pRF params on mrV >> smoothed pRF params on FS
     
     mprf_pRF_sm_FS_BS(subjid) % smoothed pRF params on FS >>  smoothed pRF params on BS
     
     % smoothed prf parameters + ROIs on BS (pial) surface saved 
     
%% Stimulus referred forward model

   % Predicted MEG response on BS  
     % input - pRF parameters in BS
     %       - MEG stimuli 
     % output - predicted MEG time series on BS
     
     %************<code>***************
     
     [prf,model] = mprf_MEG_pred_BS(subjid);
     
     % predicted response to MEG stimulus at every BS surface saved
     
   % Weighting response with gain matrix
     % input - predicted MEG time series on BS
     %       - gain matrix
     % output - predicted MEG time series on MEG sensor space
     
     %************<code>***************
     
     meg_resp = mprf_MEG_pred_Msen(subjid); 
     
     
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
   
   
  