%% mprfPreprocessMEGRetinotopyData

% This script is the main analysis to preprocess the MEG dataset.

% EK: Questions:
% - start preprocessing after which step ?. I.e. which data to load and work with?
% - In case of multiple sqd files, does data have to be combined before
% running this preprocessing script?

% addpath(genpath('~/matlab/git/toolboxes/meg_utils'))
% addpath(genpath('~/matlab/git/denoiseproject'))

% Settings: 
preproc = 'full'; % Full = do all steps,

save_interim_files = true; % Save data at every intermediate step?

% Get the time of this run:
cur_time = datestr(now);
cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';

project_dir = mprfRootPath;

cur_dir = pwd;
cd(project_dir);
% In case of synthetic data, the selected file should have an associated
% timing file, check for it and if it exist, load it and skip timing
% processing below
% [raw_file, source_dir] = uigetfile('*','Please select file to preprocess');
d = dir('/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wlsubj030/raw/R0942_RetCombined.sqd');
raw_file = d.name;
source_dir = d.folder;
put_timing_file = ['timing_info_' raw_file(5:end-4) '.mat'];

if exist(fullfile(source_dir, put_timing_file),'file')
    load(fullfile(source_dir, put_timing_file));
end
cd(cur_dir);

%source_dir = fullfile(data_dir,'Raw');
dest_dir = fullfile(source_dir,['Preprocessing_run_' cur_time]);

cd(source_dir);
mkdir(dest_dir);

%raw_file = 'R1151_Retinotopy_04.06.17.sqd';                             % Raw file

% param_dir = uigetdir(source_dir,'Select directory with parameter files');               % Parameter files
% stim_dir = uigetdir(param_dir(1:find(param_dir == '/',1,'last')),...
%     'Select directory with stimulus files');                % Stimulus files

param_dir = fullfile(source_dir, 'R0942_MegRet_9.11.18', 'behavior');
stim_dir = fullfile(source_dir, 'R0942_MegRet_9.11.18', 'stimFiles');


% Output file names:
raw_mat_data_file = 'ft_read_data_channels.mat';            % Raw mat file, read with ft_read_data, data channels only
raw_mat_tr_file = 'ft_read_trigger_channels.mat';           % Raw mat file, read with ft_read_data, trigger channels only
raw_mat_ph_file = 'ft_read_diode_channels.mat';             % Raw mat file, read with ft_read_data, photodiode only

timing_file = 'timing_info.mat';

epoched_raw_file = 'epoched_data.mat';                              % Epoched raw file
epoched_hp_filtered_file = 'epoched_data_hp.mat';                   % Epoched filtered file
epoched_hp_filtered_preproc_file = 'epoched_data_hp_preproc.mat'; %
epoched_hp_filt_preproc_denoised_file = 'epoched_data_hp_preproc_denoised.mat';




if strcmpi(preproc,'full')
    
    % Channel identity:
    trig_chan = 161 : 168;
    diode_chan = 192;
    data_chan = 1:157;
    raw_meg_file = fullfile(source_dir,raw_file);
    
    %% Triggers
    % 1. Read the triggers from the data file. This uses subfunctions that want
    % to read directly from the .sqd file. In: raw meg file, the directory
    % where all the parameter files are stored, index of the trigger and diode
    % channels. Out: timing variable that has all the trigger information you
    % need for epoching
    
    if exist('timing','var') && exist('triggers','var')
        
    else
        
        fprintf('Getting timing info...\n')
        timing = mprfGetTriggers(raw_meg_file, param_dir, trig_chan, diode_chan);
        fprintf('Done.\n');
        
        % Triggers: index in MEG signal were triggers were received and the
        % corresponding value of the trigger. Note that one point in the MEG signal
        % does not exactly correspond to one millisecond. It is 31 nanoseconds
        % longer than one millisecond. Therefore, the MEG signal seems to run
        % ahead of the stimulus/response timing, which is measured in CPU time.
        
        %TODO: make sure to deal with different outputs of mprfGetTriggers,
        %It uses many methods, some work better for one set than others.
        
        try
            triggers = [timing.trigger.trigger2flip(:,1) timing.trigger.idx(:,1)];
            
        catch
            
            triggers = [timing.trigger.channel_inds timing.trigger.idx(:)];
            
        end
        
        fprintf('Saving timing info...\n')
        save(fullfile(dest_dir,timing_file),'timing','triggers');
        fprintf('Done.\n');
    end
    
    %% Load the raw data
    % Steps:
    % 2. Read data using ft_read_data. Takes relatively long. This is the file
    % we are mainly working from, so store it separately.
    % In: raw_meg_file (.sqd) out: three files 1). All data channels (1:157)
    % .mat 2). Trigger channels (160 : 167) and 3). photodiode channel (191).
    % The reason for this is that files can get really big, so try to avoid
    % loading the entire data file when you are only interested in the
    % photodiode or trigger channels.
    
    
    
    fprintf('Reading raw data from .sqd file...\n')
    [pathstr, name, ext] = fileparts(raw_file);
    meg_raw_data = meg_load_sqd_data(source_dir, '*_Ret_*');
%     meg_raw_data = ft_read_data(raw_meg_file);
    fprintf('Done.\n');
    
    
    if save_interim_files
        
        fprintf('Saving data file...\n')
        data_channel_data = meg_raw_data(data_chan,:);
        save(fullfile(dest_dir,raw_mat_data_file), 'data_channel_data','-v7.3')
        fprintf('Done.\n');
        
        fprintf('Saving trigger file...\n')
        trigger_channel_data = meg_raw_data(trig_chan,:);
        save(fullfile(dest_dir,raw_mat_tr_file), 'trigger_channel_data','-v7.3')
        fprintf('Done.\n');
        
        fprintf('Saving photo diode file...\n')
        diode_channel_data = meg_raw_data(diode_chan,:);
        save(fullfile(dest_dir,raw_mat_ph_file), 'diode_channel_data','-v7.3')
        fprintf('Done.\n');
        
    else
        warning('Not saving intermediate data files')
        data_channel_data = meg_raw_data(data_chan,:);
        
    end
    
    
end

%% Epoching
% 3. Epoch the MEG data. Currently, epoching is done by just cutting up the
% MEG data in chunks with the length of each epoch equal to the time of the
% corresponding period. It returns a ntime_per_period by n_periods by
% n_repeats by n_channels data file. The data is padded with nans over the
% first dimension to account for differences in epoch length. These
% differences should only come about by inaccuracies of the frame timing
% during the MEG experiment.
% In: raw data, the raw data file produced by reading the .sqd file with
% ft_read_data. Can be any of the three outputs produced here, but only
% tested with the data channels

if strcmpi(preproc,'epoch') || strcmpi(preproc,'full')
    if ~exist('data_channel_data','var') || isempty(data_channel_data)
        [fname, fpath] = uigetfile('*','Please select Data channel data file');
        
        load(fullfile(fname, fpatch));
        
        if  ~exist('data_channel_data','var') || isempty(data_channel_data)
            error('Could not load data file');
            
        end
        
    end
    
    
    if ~exist('triggers','var') || isempty(triggers)
        [fname, fpath] = uigetfile('*','Please select Trigger file');
        
        load(fullefile(fname, fpatch));
        
        if  ~exist('triggers','var') || isempty(triggers)
            error('Could not load trigger file');
            
        end
        
    end
    
    epoched_data = mprfEpochMEGData(data_channel_data, triggers, [150 1100]);
    clear data_channel_data
    
    
    if save_interim_files
        fprintf('Saving epoched data...\n')
        save(fullfile(dest_dir,epoched_raw_file),'epoched_data','-v7.3');
        fprintf('Done.\n')
    else
        warning('Not saving intermediate data files')
        
    end
    
    if strcmpi(preproc, 'epoch')
        preproc = 'full';
    end
    
    
end

%% Preprocessing, high pass filtering:
% 4. High pass filtering to remove slow drift components from the MEG
% signal. This function assumes an epoched data file and concatenates the
% time series of every repeat before filtering. Thus, it does not work with
% the individual MEG epochs

if strcmpi(preproc,'filter') || strcmpi(preproc,'full')
    if ~exist('epoched_data','var') || isempty(epoched_data)
        [fname, fpath] = uigetfile('*','Please select epoched data file');
        
        load(fullfile(fpath, fname));
        
        if  ~exist('epoched_data','var') || isempty(epoched_data)
            error('Could not load data file');
            
        end
        
    end
    
    epoched_data.data = mprfFilterDataHighPass(epoched_data.data);
    
    epoched_filtered_data = epoched_data;
    
    clear epoched_data
    
    if save_interim_files
        fprintf('Saving filtered data...\n')
        save(fullfile(dest_dir,epoched_hp_filtered_file),'epoched_filtered_data','-v7.3');
        fprintf('Done.\n')
        
    else
        warning('Not saving intermediate data files')
        
    end
    
    if strcmpi(preproc, 'filter')
        preproc = 'full';
    end    
    
    
end


%% Preprocessing, reject bad epochs/channels without blink and first period epochs


if strcmpi(preproc,'preproc') || strcmpi(preproc,'full')
    
    if ~exist('epoched_filtered_data','var') || isempty(epoched_filtered_data)
        [fname, fpath] = uigetfile('*','Please select epoched filtered data file');
        
        load(fullfile(fpath, fname));
        
        if  ~exist('epoched_filtered_data','var') || isempty(epoched_filtered_data)
            error('Could not load data file');
            
        end
        
    end
    
    
    
    epoched_filtered_data = mprfPreprocessWrapper(epoched_filtered_data);
    
    epoched_filtered_preproc_data = epoched_filtered_data;
    clear epoched_filtered_data
    
    if save_interim_files
        fprintf('Saving preprocessed data...\n')
        save(fullfile(dest_dir,epoched_hp_filtered_preproc_file),'epoched_filtered_preproc_data','-v7.3');
        fprintf('Done.\n')
        
    else
        
        
    end
    
    if strcmpi(preproc, 'preproc')
        preproc = 'full';
    end
    
    
    
end

    %% Denoise the data:
    %data      : time series [channel x time samples x epoch]
    %design    : design matrix [epoch x nconds]
    

if strcmpi(preproc,'denoise') || strcmpi(preproc,'full')

    if ~exist('epoched_filtered_preproc_data','var') || isempty(epoched_filtered_preproc_data)
        [fname, fpath] = uigetfile('*','Please select epoched filtered data file');
        
        load(fullfile(fpath, fname));
        
        if  ~exist('epoched_filtered_preproc_data','var') || isempty(epoched_filtered_preproc_data)
            error('Could not load data file');
            
        end
        
    end
    
    [epoched_filtered_preproc_data, results, evalout,denoised_spec] = mprfDenoiseWrapper(epoched_filtered_preproc_data, stim_dir);
    
    epoched_filtered_preproc_denoised_data = epoched_filtered_preproc_data;
    
    clear epoched_filtered_preproc_data
    
    fprintf('Saving denoised data...\n')
    save(fullfile(dest_dir,epoched_hp_filt_preproc_denoised_file),'epoched_filtered_preproc_denoised_data','-v7.3');
    save(fullfile(dest_dir,'denoise_results'),'results','evalout','denoised_spec','-v7.3');
    
    fprintf('Done.\n')
    

        
        
       
end

return

% 
% 
% %% Preprocessing, high pass filtering:
% % 4. High pass filtering to remove slow drift components from the MEG
% % signal. This function assumes an epoched data file and concatenates the
% % time series of every repeat before filtering. Thus, it does not work with
% % the individual MEG epochs
% 
% % Load raw epoched data:
% 
% load(fullfile(data_dir,epoched_raw_file));
% 
% epoch.data = mprfFilterDataHighPass(epoch.data);
% 
% cur_time = datestr(now);
% cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
% 
% save(['epoched_data_hp_' cur_time],'epoch','-v7.3');
% 
% 
% 
% %% Preprocessing, reject bad epochs/channels without blink and first period epochs
% load(fullfile(data_dir,epoched_hp_filtered_file));
% 
% first_pos_idx = abs([0; diff(epoch.idx)]) > 0 & abs([0; diff(epoch.idx)]) < 10;
% blink_idx = epoch.idx == 20;
% 
% remove_epochs = first_pos_idx | blink_idx;
% 
% data_reshaped = epoch.data;
% data_reshaped = data_reshaped(:,~remove_epochs,:,:);
% 
% orig_dim = size(epoch.data);
% data_dim = size(data_reshaped);
% data_reshaped = reshape(data_reshaped,data_dim(1), [], data_dim(4)); % Assuming reshape does what I expect it to do....
% 
% varThreshold        = [0.05 20];
% badChannelThreshold = 0.2;
% badEpochThreshold   = 0.2;
% use3Channels        = false;
% verbose             = true;
% 
% [data_reshaped, badChannels, badEpochs]  = dfdPreprocessData(data_reshaped, ...
%     varThreshold, badChannelThreshold, badEpochThreshold,verbose);
% 
% data_reshaped(:,badEpochs,:) = nan;
% data_reshaped(:,:,badChannels) = nan;
% 
% 
% orig_data = epoch.data;
% epoch.data = nan(orig_dim);
% epoch.data(:,~remove_epochs,:,:) = reshape(data_reshaped,data_dim);
% 
% epoch.preproc.badChannels = badChannels;
% epoch.preproc.badEpochs = badEpochs;
% epoch.preproc.badChannel_thr = badChannelThreshold;
% epoch.preproc.badEpoch_thr = badEpochThreshold;
% epoch.preproc.var_thr = varThreshold;
% epoch.preproc.first_pos_idx = first_pos_idx;
% epoch.preproc.blink_idx = blink_idx;
% 
% cur_time = datestr(now);
% cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
% 
% save(['epoched_data_hp_preproc_' cur_time],'epoch','orig_data','-v7.3');
% 
% 
% 
% 
% %% Denoise the data:
% %data      : time series [channel x time samples x epoch]
% %design    : design matrix [epoch x nconds]
% clear epoch
% load(fullfile(data_dir,epoched_hp_filtered_preproc_file),'epoch');
% 
% 
% skip_n_samples = 150;
% keep_n_samples = 1100;
% fs_range = [9.9 10.1];
% 
% epoch_idx = epoch.idx;
% epoch_idx = epoch_idx(:,ones(1,19));
% epoch_idx = epoch_idx(:);
% 
% first_idx = epoch.preproc.first_pos_idx;
% blink_idx = epoch.preproc.blink_idx;
% 
% first_idx = first_idx(:,ones(1,19));
% blink_idx = blink_idx(:,ones(1,19));
% 
% first_idx = first_idx(:);
% blink_idx = blink_idx(:);
% 
% epoch_idx = epoch_idx(~first_idx &~blink_idx);
% epoch_idx = epoch_idx(~epoch.preproc.badEpochs);
% 
% data_preproc = epoch.data;
% data_preproc = reshape(data_preproc,size(data_preproc,1), [], size(data_preproc,4));
% data_preproc = permute(data_preproc, [3 1 2]);
% data_preproc = data_preproc(:,skip_n_samples+1:keep_n_samples+skip_n_samples,:);
% data_preproc = data_preproc(:,:,~first_idx & ~blink_idx);
% data_preproc = data_preproc(~epoch.preproc.badChannels,:,~epoch.preproc.badEpochs);
% 
% design = epoch_idx ~=10;
% 
% all_cond = 1:8;
% base_line = 10;
% 
% presented_conds = all_cond(ismember(all_cond,unique(epoch_idx)));
% 
% 
% stimfiles = dir(fullfile(stim_dir,'*.mat'));
% stim = load(fullfile(stim_dir, stimfiles(1).name));
% 
% stim_range = double([min(stim.stimulus.images(:)) max(stim.stimulus.images(:))]);
% bk = double(mode(mode(stim.stimulus.images(:,:,1))));
% new_period = stim.stimulus.seq(stim.stimulus.trigSeq > 0);
% 
% im_out_size = [101 101];
% im_seq = zeros([prod(im_out_size), length(new_period)],'uint8');
% 
% 
% for n = 1:length(new_period)
%     tmp_im = double(stim.stimulus.images(:,:,new_period(n)));
%     tmp_im = imresize(tmp_im,im_out_size,'nearest');
%     im_seq(:,n) = uint8(ceil(abs(tmp_im(:) - bk) ./ max(abs(stim_range - bk))));
%     
% end
% 
% design_02 = zeros(size(epoch.idx));
% n_cond = 0;
% for n = 1:length(design_02)
%     
%     if any(ismember(presented_conds, epoch.idx(n))) ...
%             && design_02(n) == 0;
%         n_cond = n_cond + 1;
%         
%         design_02(n) = n_cond;
%         
%         
%         
%         for nn = n+1:length(design_02)
%             
%             if isequal(im_seq(:,n),im_seq(:,nn)) && design_02(nn) == 0;
%                 design_02(nn) = n_cond;
%                 
%                 
%             end
%         end
%     end
%     
%     
% end
% 
% design_02 = design_02(:,ones(1,19));
% design_02 = design_02(:);
% 
% design_02 = design_02(~first_idx &~blink_idx);
% design_02 = design_02(~epoch.preproc.badEpochs);
% 
% n_conds = max(design_02);
% design_03 = zeros(size(design_02,1),n_conds);
% cond_idx = sub2ind(size(design_03),find(design_02),design_02(design_02>0));
% design_03(cond_idx) = 1;
% 
% 
% opt.verbose             = true;
% opt.use3Channels        = false;
% opt.removeFirstEpoch    = false;
% opt.removeMsEpochs      = false;
% opt.pcchoose          = 1.05;   % Take 10 PCs
% opt.npcs2try          = 10;
% 
% evokedfun = @(x)mprfDenoiseEvalFun(x,fs_range,1000);
% 
% [results,evalout,denoised_spec, denoised_data] = denoisedata(design_03,data_preproc,evokedfun,evokedfun,opt);
% 
% figure;
% plot(results.finalmodel.r2,results.origmodel.r2,'k.');
% hold on;
% plot([0 10],[0 10],'r-');
% xlabel('Final');
% ylabel('orig');
% 
% orig_dim = [size(denoised_data{1},2) 2660 157];
% 
% first_pos_idx = epoch.preproc.first_pos_idx;
% blink_idx = epoch.preproc.blink_idx;
% 
% first_pos_idx = first_pos_idx(:,ones(1,19));
% blink_idx = blink_idx(:,ones(1,19));
% 
% first_pos_idx = first_pos_idx(:);
% blink_idx = blink_idx(:);
% 
% all_idx = first_pos_idx | blink_idx;
% 
% good_epoch_idx = false(size(all_idx));
% good_epoch_idx(~all_idx) = ~epoch.preproc.badEpochs;
% 
% tmp_data = nan(orig_dim);
% tmp_data(:,good_epoch_idx,~epoch.preproc.badChannels) = permute(denoised_data{1},[2 3 1]);
% 
% epoch.data = reshape(tmp_data,[1100 140 19 157]);
% 
% cur_time = datestr(now);
% cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
% 
% save(['epoched_hp_filtered_preprocessed_denoised_data_' cur_time],'epoch','-v7.3');
% save(['denoise_output_' cur_time],'results','denoised_spec','evalout');
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
















