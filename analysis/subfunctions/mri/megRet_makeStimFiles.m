function [] = makeStimFiles(subject, sessionDir, sessionName, bidsSession)

% Create scan_params and scan_images matfiles to run retinotopy model

% Where to find stimuli params?
d1 = dir(fullfile(sessionDir, 'Stimuli','stimFiles', '*stimulus*'));
d2 = dir(fullfile(sessionDir, 'Stimuli','paramsFiles', '*.mat'));

% Where to save stimuli params?
stimFolder = fullfile(sessionDir, 'Stimuli');
if ~exist(stimFolder, 'dir'); mkdir(stimFolder); end;


% Just pick the first run:
id = 1;

% create params struct
load(fullfile(d2(id).folder, d2(id).name));

% params.analysis.fieldSize   = 11.18;
% params.analysis.numberStimulusGridPoints = 11;
% params.analysis.sampleRate  = 1.0538;
% params.prescanDuration      = 8; % TRs

    
% im_fName = fullfile(d1(id).folder, d1(id).name);
% pars_fName = fullfile(d2(id).folder, d2(id).name);
% 
% params.stim(id).imFile          = im;
% params.stim(id).paramsFile      = pars;
% params.stim(id).framePeriod     = 1.5; % seconds
% params.stim(id).nFrames         = 244; % TRs
% params.stim(id).imFilter        = 'binary';
%     
% params = makeStimFromScan(params,id);


% add additional parameters
im = load(fullfile(d1(id).folder, d1(id).name),'stimulus');
params.seq          =  im.stimulus.seq;
params.seqtiming    =  im.stimulus.seqtiming;
params.ncycles      =  1;
params.numCycles    =  params.ncycles;
if strcmp(subject, 'wlsubj030')
    params.prescanDuration = 8*1.5; % Sec
    params.scanDuration =  244-8; % TRs
elseif strcmp(subject, 'wlsubj068')
    params.prescanDuration = 0; %12*1.5; % Sec
    params.scanDuration =  244; %244-12; % TRs
end

params.period       =  params.scanDuration;
images              = im.stimulus.images;

save(fullfile(stimFolder,'scan_images.mat'), 'images')

clear im images

save(fullfile(stimFolder,'scan_params.mat'));
    



% id = id: scan number (integer);

%
%       params: a struct containing the following fields
%           (normally generated from GUI or stored in dataTYPES)
%
%          Required
%               params.stim(id).imFile*  (filename)
%               params.stim(id).paramsFile** (filename)
%               params.analysis.fieldSize (degrees, radius)
%               params.analysis.numberStimulusGridPoints (n points, radius)
%               params.analysis.sampleRate (degrees per point)
%
%          Optional
%               params.framePeriod (length of TR in s)
%               params.stim(id).nFrames (n TRs in scan)
%               params.prescanDuration (in TRs, not s)
%               params.stim(id).imFilter (default = 'binary');
%                   (see .../retinotopyModel/FilterDefinitions/ for 
%                       other filters)
%
%
%       *imFile:  imfile wil be loaded into struct 'I'
%           Required
%               I.images: an n x m x k matrix. n and m are the image size,
%                           k is the number of unique images.
%                TODO: allow RGB image matrices 
%
%       **paramsFile: paramsfile wil be loaded into struct 'P'.
%           Required
%               P.stimulus.seq: a vector indexing the image matrix I.images
%               P.stimulus.seqTiming: a vector of image onset times (in s)
%
%           Optional
%               P.params.display (screen calibration information)
%           Optional (but nec if not included in input params)
%               P.params.framePeriod (length of TR in s)
%               P.params.numImages (n TRs in scan)
%               P.params.prescanDuration (in s, not TRs)








% OBSOLETE:
% 
% % Check size images:
% disp(size(stimulus.images))
% 
% % Allocate space for images array
% images = zeros(size(stimulus.images,1),size(stimulus.images,2),length(unique(stimulus.seq)));
% 
% % unique stim
% u_stim = unique(stimulus.seq);
% 
% % Loop over sequence of images 
% for ii = 1:length(u_stim)
%     
%     images(:,:,ii) =  stimulus.images(:,:,stimulus.seq(u_stim(ii)));
%     
% end
% 
% 
% save(fullfile(stim_folder,'scan_images.mat'), 'images', '-v7.3');
% 
% % for visualization:
% figure; for jj = 1:size(images,3); imagesc(images(:,:,jj)); colormap gray; drawnow; title(jj); pause(0.1); end;
% 
% %% Save scan params
% 
% clear all;
% 
% % Where does retinotopy param file live:
% load('/Volumes/server/Projects/MEG/Retinotopy/Data/fMRI/wlsubj030/Stimuli/vistadisp_files/vistadisp_output_matrices/SL1_20180910T103404.mat')
% 
% prescanDuration = 8*1.5; % 8 volumes pre blank, with a TR of 1.5s;
% nCycles         = 1;
% 
% 
% session_path = '/Volumes/server/Projects/MEG/Retinotopy/Data/fMRI/wlsubj030';
% stim_folder = fullfile(session_path, 'Stimuli');
% 
% save(fullfile(stim_folder,'scan_params.mat'));
