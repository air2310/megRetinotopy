function mprfSession_preprocess_meg_data(do_syn, syn_fname_in)

if ~exist('do_syn','var') || isempty(do_syn)
    do_syn = false;
end

if do_syn
    if ~exist('syn_fname_in','var') || isempty(syn_fname_in)
        error('Need synthetic data file name')
    end
    
    if ~exist(syn_fname_in,'file')
        error('Synthetic data file does not exist')
    end
    
end


global mprfSESSION

if all(size(mprfSESSION) == [0 0])
    load('mprfSESSION.mat')
    
end

main_dir = mprf__get_directory('main_dir');
meg_ddir = fullfile(main_dir, mprf__get_directory('meg_data'));
syn_ddir = fullfile(main_dir, mprf__get_directory('syn_data'));
stim_dir = fullfile(main_dir, mprf__get_directory('meg_imported_stim'));

has_meg_file = false;
has_syn_file = false;
has_stim_file = false;

meg_files = {};
syn_files = {};
stim_files = {};

if do_syn
    
    syn_files = {syn_fname_in};
    has_syn_file = true;
else
    if exist(meg_ddir,'dir')
        in_meg_dir = dir(fullfile(meg_ddir, '*.sqd'));
        if ~isempty(in_meg_dir)
            has_meg_file = true;
            meg_files = {in_meg_dir.name};
            
        end
    end
    
    if exist(syn_ddir,'dir')
        in_syn_dir = dir(fullfile(syn_ddir, '*.sqd'));
        if ~isempty(in_syn_dir)
            has_syn_file = true;
            syn_files = {in_syn_dir.name};
            
            
        end
    end
    
end

if exist(stim_dir,'dir')
    in_stim_dir = dir(fullfile(stim_dir, '*.mat'));
    if ~isempty(in_stim_dir)
        has_stim_file = true;
        stim_files = {in_stim_dir.name};
        
    end
end

if (has_meg_file || has_syn_file) && has_stim_file
    
else
    fprintf('Error, no data and/or stimulus found\n')
    return
end


params = mprf__preprocess_meg_gui(meg_files, syn_files, stim_files, do_syn);


save_interim_files = params.preproc.do.save; % Save data at every intermediate step?
% Get the time of this run:
cur_time = datestr(now);
cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';

% Main paths:
if ~exist('denoisedata','file')
    tbUse({'retMEG','knk'})
    add_dirs('denoise')
    
end


if strcmpi(params.preproc.data.type,'meg_data');
    raw_meg_file = fullfile(main_dir, mprf__get_directory('meg_data'),params.preproc.data.file);
    [~,fname] = fileparts(raw_meg_file);
    fname(fname == ' ' | fname == ':' | fname == '-'|fname == '.') = '_';
    dest_dir = fullfile(main_dir, mprf__get_directory('meg_preproc'),['pp_run_' fname '_on_' cur_time]);
    
    mkdir(dest_dir)
    
elseif strcmpi(params.preproc.data.type,'syn_data');
    raw_meg_file = fullfile(main_dir, mprf__get_directory('syn_data'),params.preproc.data.file);
    [fpath,fname] = fileparts(raw_meg_file);
    fname(fname == ' ' | fname == ':' | fname == '-'|fname == '.') = '_';
    dest_dir = fullfile(main_dir, mprf__get_directory('syn_preproc'),['pp_run_' fname '_on_' cur_time]);
    
    tmp = strfind(params.preproc.data.file,'syn_data_');
    if any(tmp)
        tmp2 = params.preproc.data.file(tmp+length('syn_data_'):end);
        [~, fname] = fileparts(tmp2);
        trig_info_name_guess = fullfile(fpath,['trig_info_' fname '.mat']);
        
        if exist(trig_info_name_guess,'file')
            load(trig_info_name_guess)
            fprintf('Loading %s as timing file\n',trig_info_name_guess);
        end
    end
    
    mkdir(dest_dir)
    
end

param_dir = fullfile(main_dir,mprf__get_directory('meg_stim_params'));
stim_dir = fullfile(main_dir, mprf__get_directory('meg_stim'));


% Output file names:
raw_mat_data_file = 'ft_read_data_channels.mat';            % Raw mat file, read with ft_read_data, data channels only
raw_mat_tr_file = 'ft_read_trigger_channels.mat';           % Raw mat file, read with ft_read_data, trigger channels only
raw_mat_ph_file = 'ft_read_diode_channels.mat';             % Raw mat file, read with ft_read_data, photodiode only

timing_file = 'timing_info.mat';

epoched_raw_file = 'epoched_data.mat';                              % Epoched raw file
epoched_hp_filtered_file = 'epoched_data_hp.mat';                   % Epoched filtered file
epoched_hp_filtered_preproc_file = 'epoched_data_hp_preproc.mat'; %
epoched_hp_filt_preproc_denoised_file = 'epoched_data_hp_preproc_denoised.mat';


% Channel identity:
trig_chan = params.channels.triggers;
diode_chan = params.channels.diode;
data_chan = params.channels.data;

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
        triggers = [timing.trigger.trigger2flip(:,1) timing.trigger.idx(:)];
        
    catch
        
        triggers = [timing.trigger.channel(:,1) timing.trigger.idx(:)];
        
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
meg_raw_data = ft_read_data(raw_meg_file);
fprintf('Done.\n');


if save_interim_files
    
    fprintf('Saving data file...\n')
    data_channel_data = meg_raw_data(data_chan+1,:);
    save(fullfile(dest_dir,raw_mat_data_file), 'data_channel_data','-v7.3')
    fprintf('Done.\n');
    
    fprintf('Saving trigger file...\n')
    trigger_channel_data = meg_raw_data(trig_chan+1,:);
    save(fullfile(dest_dir,raw_mat_tr_file), 'trigger_channel_data','-v7.3')
    fprintf('Done.\n');
    
    fprintf('Saving photo diode file...\n')
    diode_channel_data = meg_raw_data(diode_chan+1,:);
    save(fullfile(dest_dir,raw_mat_ph_file), 'diode_channel_data','-v7.3')
    fprintf('Done.\n');
    
else
    warning('Not saving intermediate data files')
    data_channel_data = meg_raw_data(data_chan+1,:);
    
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

if ~exist('data_channel_data','var') || isempty(data_channel_data)
    [fname, fpath] = uigetfile('*','Please select Data channel data file');
    
    load(fullefile(fname, fpatch));
    
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

epoched_data = mprfEpochMEGData(data_channel_data, triggers, [params.preproc.params.epoch.low params.preproc.params.epoch.include]);
clear data_channel_data


if save_interim_files
    fprintf('Saving epoched data...\n')
    save(fullfile(dest_dir,epoched_raw_file),'epoched_data','-v7.3');
    fprintf('Done.\n')
else
    warning('Not saving intermediate data files')
    
end


%% FILTER THE DATA

if params.preproc.do.filt
    if ~exist('epoched_data','var') || isempty(epoched_data)
        [fname, fpath] = uigetfile('*','Please select epoched data file');
        
        load(fullfile(fpath, fname));
        
        if  ~exist('epoched_data','var') || isempty(epoched_data)
            error('Could not load data file');
            
        end
        
    end
    
    epoched_data.data = mprfFilterDataHighPass(epoched_data.data);
    
    if save_interim_files
        fprintf('Saving filtered data...\n')
        save(fullfile(dest_dir,epoched_hp_filtered_file),'epoched_data','-v7.3');
        fprintf('Done.\n')
        
    else
        warning('Not saving intermediate data files')
        
    end
    
    
    
end


%% Preprocessing, reject bad epochs/channels without blink and first period epochs


if params.preproc.do.preproc
    
    if ~exist('epoched_data','var') || isempty(epoched_data)
        [fname, fpath] = uigetfile('*','Please select epoched data file to process');
        
        load(fullfile(fpath, fname));
        
        if  ~exist('epoched_filtered_data','var') || isempty(epoched_data)
            error('Could not load data file');
            
        end
        
    end
    
    
    
    epoched_data = mprfPreprocessWrapper(epoched_data);
    
    if save_interim_files
        fprintf('Saving preprocessed data...\n')
        save(fullfile(dest_dir,epoched_hp_filtered_preproc_file),'epoched_data','-v7.3');
        fprintf('Done.\n')
        
    else
        
        
    end
    
end


%% Denoise the data:
%data      : time series [channel x time samples x epoch]
%design    : design matrix [epoch x nconds]

if params.preproc.do.denoise
    
    if ~exist('epoched_data','var') || isempty(epoched_data)
        [fname, fpath] = uigetfile('*','Please select epoched data to process');
        
        load(fullfile(fpath, fname));
        
        if  ~exist('epoched_data','var') || isempty(epoched_data)
            error('Could not load data file');
            
        end
        
    end
    
    [epoched_data, results, evalout,denoised_spec] = mprfDenoiseWrapper(epoched_data, stim_dir);
    
    fprintf('Saving denoised data...\n')
    
    if save_interim_files
        save(fullfile(dest_dir,epoched_hp_filt_preproc_denoised_file),'epoched_data','-v7.3');
        
    end
    save(fullfile(dest_dir,'denoise_results'),'results','evalout','denoised_spec','-v7.3');
    fprintf('Done.\n')
    
end

if ~save_interim_files
    out_name = 'epoched_data';
    field_names = fieldnames(params.preproc.do);
    
    for n  = 1:length(field_names)
        cur_field = field_names{n};
        cur_val = params.preproc.do.(cur_field);
        switch lower(cur_field)
            
            case 'filt'
                if cur_val
                   out_name = [out_name '_hp'];
                    
                end
                
            case 'preproc'
                 if cur_val
                   out_name = [out_name '_preproc '];
                    
                end
                
                
                
            case 'denoise'
                if cur_val
                    out_name = [out_name '_denoised '];
                    
                end
                
                
        end
    end
    save(fullfile(dest_dir,[out_name, '.mat']),'epoched_data','-v7.3');
    
    

end



