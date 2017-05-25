
% Keep track of directories:
cur_dir = pwd;

% Look for these triggers:
cond_triggers = [4 5 6 10 12 15 20 30];

% Trigger channels and diode channels in the MEG data file. Zero based.
trig_chan = 160 : 167;
diode_chan = 191;

% Directory where the raw data is located:
raw_dir = '/Users/bpklein/Documents/Server_offline_folder/Data/MEG/wl_subj040/Ret_check_01/Raw';

% Name of data file to load:
fname = 'R1151_RetCheckerboard1_3.17.17.sqd';

% All files in the parameter directory:
param_dir = '/Users/bpklein/Documents/Server_offline_folder/Data/MEG/wl_subj040/R11551_Retinotopy_3.17.17';
par_files = dir(fullfile(param_dir,'*.mat'));

% Load the data:
[data, inf] = sqdread(fullfile(raw_dir, fname),'Channels',[trig_chan diode_chan]);

% Triggers:
tdata = data(:,1:8);

% Work out the trigger types and their timing:
tmp_triggers = meg_fix_triggers(tdata);
tmp = find(tmp_triggers);
triggers = [tmp(:) tmp_triggers(tmp)];


% Photo diode:
phdata = data(:,9);

clear data tdata;



% How many times is each condition repeated:
nrepeats = sum(bsxfun(@rdivide,bsxfun(@times, triggers(:,2), cond_triggers),cond_triggers.^2) == 1);

% Throw an error for now when the amount of repeats is not the same for
% every condition:
if length(unique(nrepeats)) == 1;
    nrepeats = nrepeats(1);
else
    error('Not all conditions are repeated the same amount of times');
end

% Get the order in which the conditions are presented for every repeat.
% This also throws an error when the amount of repeats is not the same for
% every condtion:
cond_order = reshape(triggers(ismember(triggers(:,2),cond_triggers),2),length(cond_triggers),nrepeats);
file_idx = nan(1,nrepeats);
% first_flip_times = nan(1,nrepeats);
% last_flip_times = nan(1,nrepeats);
%
% flip_idx = cell(1,nrepeats);
% stim_trig_idx = cell(1,nrepeats);
% trig_flip_idx = cell(1,nrepeats);
% trig_seq_time = cell(1,nrepeats);
% trig_flip_time = cell(1,nrepeats);

% Loop over the parameter files and find matching condition orders:
for n = 1:length(par_files)
    load(fullfile(param_dir,par_files(n).name),'stimulus');
    [aa, bb] = unique(stimulus.trigSeq);
    aa2 = aa(ismember(aa, cond_triggers));
    bb2 = bb(ismember(aa, cond_triggers));
    
    [~,idx] = sort(bb2);
    
    cur_cond_order = aa2(idx);
    cond_match = find(mean(bsxfun(@eq, cond_order, cur_cond_order)) == 1);
    
    
    if cond_match
        if length(cond_match) > 1 % In this case, there are more matches, in the end, keep both files and work the correct file out by assessing their timing
            error('Conditions are not unique, can not resolve to which repeat which file belongs')
            
        else
            file_idx(cond_match) = n;
            fprintf('Matching file found\n');
            
            
            
        end
    end
    
    
end

time0 = [];
for n = 1:length(file_idx)
    
    load(fullfile(param_dir,par_files(file_idx(n)).name),'response','stimulus','params'); % Also need stimulus
    [flip_time, stim_time, init, time0] = mprfGetStimAndFlipTime(stimulus,response, params, time0);

    
    if n == 1;
        stim_seq_time = nan(size(stim_time.seq_times,1), length(file_idx));
        stim_flip_time = nan(size(flip_time.flip_times,2), length(file_idx));
        trigger_flip_time = nan(size(flip_time.trigger_times,2), length(file_idx));
        trigger_seq_time = nan(size(stim_time.trigger_times,1), length(file_idx));
        trigger_flip_time_02 = nan(size(flip_time.trigger_times_02,1), length(file_idx));
    end
    
    stim_seq_time(:,n) = stim_time.seq_times; % Time at which a flip was requested
    stim_flip_time(:,n) = flip_time.flip_times'; % Timing of the flips, relative to start of the experiment
    trigger_flip_time(:,n) = flip_time.trigger_times'; % time at which a trigger occured
    trigger_seq_time(:,n) = stim_time.trigger_times; % time at which a trigger was requested
    trigger_flip_time_02(:,n) = flip_time.trigger_times_02; % trigger value for each flip, zeros is no trigger
    
    
    
    
end

stim_seq_time = stim_seq_time(:);
stim_flip_time = stim_flip_time(:);

delay_sequence = stim_flip_time - stim_seq_time; % Difference between sequence and flips
d_delay_sequence = [0; diff(delay_sequence)]; % Differential, the flips that were delayed

figure;
plot(delay_sequence)
ylabel('Delay (milliseconds)');


figure; 
plot(stim_flip_time(:), d_delay_sequence)

trigger_interval = diff(round(triggers(:,1)) ./ 1000);
flip_interval = diff(trigger_flip_time(:) ./ 1000);
seq_interval = diff(trigger_seq_time(:) ./ 1000);

% Amount of frames the triggers are lagging behind the flips:
trigger.trigger_flip_delay = round((flip_interval - trigger_interval) ./ (1/params.display.frameRate));

% Amount of frames the flips are lagging behind the stimulus sequence
trigger.flip_seq_delay = round((flip_interval - seq_interval) ./ (1/params.display.frameRate));

% Amount of frames the triggers are lagging behind the stimulus sequence
trigger.trigger_seq_delay = round((trigger_interval - seq_interval) ./ (1/params.display.frameRate));

accum_delay = 0;
delay_idx = nan(1,10^5);
delay_idx_start = 1;
for n = 1:length(stim_flip_time)
    
    cur_delay = stim_flip_time(n) - stim_seq_time(n);
    added_delay = cur_delay - accum_delay;
    accum_delay = accum_delay + added_delay;
    
    if added_delay > floor(1/params.display.frameRate * 1000);
        delay_start = stim_seq_time(n) + (stim_seq_time(n) - stim_seq_time(n-1))+1;
        delay_end =  delay_start + added_delay-1;
        
        delay_idx_end = delay_idx_start + added_delay-1;
        delay_idx(delay_idx_start:delay_idx_end) = delay_start : delay_end;
        
        
        delay_idx_start = delay_idx_start + added_delay;
        
    end
    
end

delay_idx = delay_idx(~isnan(delay_idx)); % indices of the time points that 
%were delayed. I.e. where the flip timing and sequence timing mismatch












