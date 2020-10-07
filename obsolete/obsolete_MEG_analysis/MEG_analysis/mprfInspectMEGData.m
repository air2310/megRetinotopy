% Goal: plot triggers, trigger sequence, stimulus sequence and photodiode
% sequence for comparison. Needed: all param files associated with the
% current MEG time course. The time at which the files are created may be
% useful. But it is probably better to only check this manually. Directly
% after the experiment organize the param files... For now, check how many
% stimulus repetitions are presented. Easy to by assessing how much repeats
% of the stimulus specific triggers there are.

% Issues: The photodiode flip times (in MEG signal time) worked out from
% the start of the init sequence in the MEG photodiode trace and the first
% flip time of this init squence, which gives the MEG time corresponding to this flip time.
% I add the delay in flip time of the diode triggers to this number (in MS)
% to get the diode flip in MEG time. This already incorporates the delay
% between the diode flip in CPU time and MEG time. We simply do not have 
% from the parameter file of the corresponding run. We simply do not know
% the time at which MEG data collection started, so aliging the stimulus
% sequence/photodiode sequence/triggers etc to the MEG signal this way is
% already delayed. 



% Keep track of directories:
cur_dir = pwd;

% Look for these triggers:
cond_triggers = [4 5 6 10 12 15 20 30];

% Trigger channels and diode channels in the MEG data file. Zero based.
trig_chan = 160 : 167;
diode_chan = 191;

% Directory where the raw data is located:
raw_dir = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wlsubj030/raw';

% Name of data file to load:
fname = 'R0942_Ret2_9.11.18.sqd';

% All files in the parameter directory:
param_dir = fullfile(raw_dir,'R0942_MegRet_9.11.18','behavior');
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

if sum(isnan(file_idx)) == 0
else
    error('Could not find a matching file for every repeat');
    
end


% Process trigger data:
% rescale to [0 1]
phdata = phdata - min(phdata(:));
phdata = phdata / max(phdata(:));

% check whether triggers are indicated by a low value or a high value
if round(mean(phdata)) == 0, trigger_is_high = true;
else                     trigger_is_high = false; end

% if triggers are indicated by a low value, then invert the signal
if ~trigger_is_high, phdata = 1 - phdata; end

% threshold to binarize from analog recording
phdata = phdata > 0.1;


% differentiate to isolate trigger onsets, and pad to preserve length
phdata_d = [diff(phdata);0];

% rectify
phdata_d(phdata_d<0) = 0;

% mark the samples as 1 when any trigger channel signaled
any_trigger      = sum(phdata_d,2) > 0;

% list the sample numbers (time ponts) when any trigger channel signalled
any_trigger_inds = find(any_trigger);


npoints_bin = 1000;

phdata_d = double(phdata_d);
bin = zeros(npoints_bin,1);
bin_vector = zeros(size(phdata_d));

trial_start_points = nan(nrepeats);
trial_end_points = nan(nrepeats);


stim_start_points = nan(nrepeats);
diode_flip_times = nan(16, nrepeats);

for n = 1:length(file_idx)
    
    load(fullfile(param_dir, par_files(file_idx(n)).name),'params','stimulus','response','time0');
%    [~,t_saved] = fileparts(par_files(file_idx(n)).name)
    
    init_seq = params.display.initstim.seq;     % The sequence used to intialize the diode, i.e when is white (1) when is black (0)
    %init_seqtiming = params.display.initstim.seqtiming; % Timing of the diode flashes,
    % They are relative to the previous flip of the previous flash
    
    init_seq_cpu_time  = stimulus.flashTimes.flip; % CPU time at which flips of flashes occured
    init_seq_cpu_flash_times = init_seq_cpu_time(logical(init_seq));

    % CPU TIME:
    % Delay last flash and frist stimulus presentation:
    time_first_flash = init_seq_cpu_flash_times(1);
    init_duration = round((init_seq_cpu_flash_times(end) - time_first_flash) * 1000);
    time_first_frame = response.flip(1);
    
    flash_stim_delay = round((time_first_frame - time_first_flash) * 1000); % ms
    diode_flip_delays = round((response.flip(logical(stimulus.diodeSeq)) - time_first_flash) * 1000);
    
    %trigger_sep_ms = diff(any_trigger_inds) / 1000;
    %sep_between_diode_triggers = reshape(trigger_sep_ms([2:2+3 25:25+3 48:48+3 71:71+3 94:94+3]),4,5); % Separation between diode triggers
    
    sep_between_flashes = [0 cumsum(round(diff(init_seq_cpu_flash_times) * 1000))];
    
    sep_ind = [sep_between_flashes]+1;
    
    
    npoints_bin = 1000;
    bin = zeros(npoints_bin,1);
    bin(sep_ind(:)) = 1;
    
    bin_vector((length(phdata_d)-npoints_bin)/2+1:(length(phdata_d)-npoints_bin)/2+npoints_bin) = bin;
    
    ft_bin =fft(bin_vector);
    ft_ph = fft(double(phdata_d));
    
    aa = ifftshift(ifft(ft_bin .* ft_ph));
%     figure; plot(aa);
    
    trial_start_points(:,n) = find(aa == max(aa));
    trial_end_points(:,n) = trial_start_points(:,n) + init_duration;
    stim_start_points(:,n) = trial_start_points(:,n) + flash_stim_delay;
    diode_flip_times(:,n) = trial_start_points(n,n) + diode_flip_delays;
    
end

trial_start = trial_start_points(logical(eye(nrepeats)));
stim_start =   stim_start_points(logical(eye(nrepeats)));

to_check = zeros(size(phdata_d));
to_check(trial_start) = 0.22;
to_check(stim_start) = .24;

trigger_check = zeros(size(phdata_d));
trigger_check(triggers(:,1)) = 0.18;

diode_trigger_check = zeros(size(phdata_d));
diode_trigger_check(diode_flip_times(:)) = 0.18;

plot_idx = 1:50000;

figure;
plot(to_check(plot_idx));
hold on;
plot(phdata_d(plot_idx) ./ 5,'r')
plot(trigger_check(plot_idx),'k')
plot(diode_trigger_check(plot_idx),'y')


diode_flip_times2 = diode_flip_times(:); % These are the times at which the flips occured of the diode sequence.
% obtained by adding the fliptimes of the triggers to the estimated start
% of the stimulus in the MEG signal (estimated by the autocorrelation above).

% Compute the time difference between the diode flip times and the times
% at which a flash was detected:
ph_trig_diff = bsxfun(@minus, diode_flip_times2, any_trigger_inds');

% Select the flips and triggers that are closest to each other:
[min_diffs, min_idx] = min(abs(ph_trig_diff),[],1);
[min_diffs2, min_idx2] = min(abs(ph_trig_diff),[],2);

idx = sub2ind(size(ph_trig_diff),min_idx,1:size(ph_trig_diff,2));
idx2 = sub2ind(size(ph_trig_diff),1:size(ph_trig_diff,1),min_idx2');

% Remove any trigger that was more than 100 ms away from a photodiode
% flip, i.e. we simply did not send a photodiode flip for this
% trigger (i.e. they are from the init sequence):
keep_idx = min_idx(min_diffs < 100);
keep_idx2 = min_idx2(min_diffs2 < 100);

% These triggers were matched by the corresponding photodiode triggers:
matched_flips = diode_flip_times2(keep_idx);
matched_triggers = any_trigger_inds(keep_idx2);

figure;
plot(matched_flips - matched_triggers,'k.-')
title('Delay between diode flip and corresponding trigger in sequence')
ylabel('Milliseconds')
xlabel('N trigger')

return
cur_dir = pwd;

% Channels are zero based:
trig_chan = 160 : 167;
diode_chan = 191;

raw_dir = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wl_subj040/Ret_check_01/Raw';
fname = 'R1151_RetCheckerboard1_3.17.17.sqd';
param_dir = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wl_subj040/R11551_Retinotopy_3.17.17';
par_files = dir(param_dir);

[tdata, inf] = sqdread(fullfile(raw_dir, fname),'Channels',trig_chan);
phdata =  sqdread(fullfile(raw_dir, fname),'Channels',diode_chan);

tmp_triggers = meg_fix_triggers(tdata);

tmp = find(tmp_triggers);
triggers = all_trigger(fullfile(raw_dir, fname),trig_chan);
triggers_02 = [tmp(:) tmp_triggers(tmp)];


% rescale to [0 1]
phdata = phdata - min(phdata(:));
phdata = phdata / max(phdata(:));

% check whether triggers are indicated by a low value or a high value
if round(mean(phdata)) == 0, trigger_is_high = true;
else                     trigger_is_high = false; end

% if triggers are indicated by a low value, then invert the signal
if ~trigger_is_high, phdata = 1 - phdata; end

% threshold to binarize from analog recording
phdata = phdata > 0.1;


% differentiate to isolate trigger onsets, and pad to preserve length
phdata_d = [diff(phdata);0];

% rectify
phdata_d(phdata_d<0) = 0;

% mark the samples as 1 when any trigger channel signaled
any_trigger      = sum(phdata_d,2) > 0;

% list the sample numbers (time ponts) when any trigger channel signalled
any_trigger_inds = find(any_trigger);

% Compute the difference between the trigger timpe points and the
% photodiode timepoints:
ph_trig_diff = bsxfun(@minus, triggers(:,1), any_trigger_inds');
ph_trig_diff2 = bsxfun(@minus, triggers_02(:,1), any_trigger_inds');

% Select the triggers that are closest to a photodiode trigger:
[min_diffs, min_idx] = min(abs(ph_trig_diff),[],1);
[min_diffs2, min_idx2] = min(abs(ph_trig_diff),[],2);

% Remove any trigger that was more than 100 ms away from a photodiode
% trigger, i.e. we simply did not send a photodiode trigger for this
% trigger:
keep_idx = min_idx(min_diffs < 100);
keep_idx2 = min_idx2(min_diffs2 < 100);

% These triggers were matched by the corresponding photodiode triggers:
matched_triggers = triggers(keep_idx,:);
matched_diodes = any_trigger_inds(keep_idx2);

% METHOD 2


% Select the triggers that are closest to a photodiode trigger:
[min_diffs, min_idx] = min(abs(ph_trig_diff2),[],1);
[min_diffs2, min_idx2] = min(abs(ph_trig_diff2),[],2);

% Remove any trigger that was more than 100 ms away from a photodiode
% trigger, i.e. we simply did not send a photodiode trigger for this
% trigger:
keep_idx = min_idx(min_diffs < 100);
keep_idx2 = min_idx2(min_diffs2 < 100);

% These triggers were matched by the corresponding photodiode triggers:
matched_triggers_02 = triggers_02(keep_idx,:);
matched_diodes_02 = any_trigger_inds(keep_idx2);



[matched_diodes matched_triggers matched_diodes_02 matched_triggers_02];



% From the retinotopy code
load(fullfile(param_dir, par_files(3).name));
[~,t_saved] = fileparts(par_files(3).name);

init_seq = params.display.initstim.seq;     % The sequence used to intialize the diode, i.e when is white (1) when is black (0)
init_seqtiming = params.display.initstim.seqtiming; % Timing of the diode flashes,
% They are relative to the previous flip of the previous flash

init_seq_cpu_time  = stimulus.flashTimes.flip; % CPU time at which flips of flashes occured
int_seq_cpu_flash_times = init_seq_cpu_time(logical(init_seq));


%trigger_sep_ms = diff(any_trigger_inds) / 1000;
%sep_between_diode_triggers = reshape(trigger_sep_ms([2:2+3 25:25+3 48:48+3 71:71+3 94:94+3]),4,5); % Separation between diode triggers

sep_between_flashes = [0 cumsum(round(diff(int_seq_cpu_flash_times) * 1000))];

sep_ind = [sep_between_flashes]+1;


npoints_bin = 1000;
bin = zeros(npoints_bin,1);
bin(sep_ind(:)) = 1;

bin_vector = [zeros((length(phdata)-npoints_bin)/2,1) ; bin ; zeros((length(phdata)-npoints_bin)/2,1)];

ft_bin =fft(bin_vector);
ft_ph = fft(double(phdata_d));

aa = ifftshift(ifft(ft_bin .* ft_ph));
figure; plot(aa);
tmp = floor(npoints/2);

find(phdata_d == 1,100,'first');
find(aa == 3);



