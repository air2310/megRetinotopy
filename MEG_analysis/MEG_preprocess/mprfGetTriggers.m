function timing = mprfGetTriggers(meg_file, param_dir, trig_chan, diode_chan, do_plot)
if ~exist('do_plot','var') || isempty(do_plot)
    do_plot = false;
end

param_files = dir(fullfile(param_dir,'*.mat')); % Make sure this is in chronological order. I.e. first file in struct corresponds to the first stimulus run


%% Reconstruct trigger sequence from stimulus/response information...

time0 = [];
for n = 1:length(param_files)
    
    load(fullfile(param_dir,param_files(n).name),'response','stimulus','params'); % Also need stimulus
    [flip_time, stim_time, init, time0] = mprfGetStimAndFlipTime(stimulus,response, params, time0);
    
    
    if n == 1;
        timing.stimulus.seqtime = nan(size(stim_time.seq_times',2), length(param_files));
        timing.stimulus.flip_time = nan(size(flip_time.flip_times,2), length(param_files));
        
        timing.trigger.flip_time = nan(size(flip_time.trigger_times,2), length(param_files));
        timing.trigger.seq_time = nan(size(stim_time.trigger_times',2), length(param_files));
        timing.trigger.flip_time_02 = nan(size(flip_time.trigger_times_02',2), length(param_files));
        
        timing.trigger.idx = nan(size(flip_time.trigger_times,2), length(param_files));
        timing.init.flip_time = nan(size(init.flip_times,2), length(param_files));
        timing.init.idx = nan(size(init.seq.seq,2), length(param_files));
        timing.init.flash_time = nan(sum(init.seq.seq), length(param_files));
        
    end
    
    timing.stimulus.seqtime(:,n) = stim_time.seq_times'; % Time at which a flip was requested
    timing.stimulus.flip_time(:,n) = flip_time.flip_times'; % Timing of the flips, relative to start of the experiment
    
    if ~isempty(flip_time.trigger_times)
        timing.trigger.flip_time(:,n) = flip_time.trigger_times'; % time at which a trigger occured
    end
    
    if ~isempty(flip_time.trigger_times)
        timing.trigger.seq_time(:,n) = stim_time.trigger_times'; % time at which a trigger was requested
    end
    
    timing.trigger.flip_time_02(:,n) = flip_time.trigger_times_02; % trigger value for each flip, zeros is no trigger
    timing.trigger.idx(:,n) = stim_time.trigger_idx';
    
    timing.init.flip_time(:,n) = init.flip_times;
    timing.init.idx(:,n) = init.seq.seq;
    timing.init.flash_time(:,n) = init.flip_times(logical(init.seq.seq));
    
    if do_plot
        figure;
        
        % desired inter-stimulus duration
        plot(diff(stimulus.seqtiming));
        
        % measured inter-stimulus duration
        hold on; plot(diff(response.flip), 'r-');
        
        ylim(median(diff(response.flip)) + [-.001 .001])
        % frames between stimuli
        frames = round(diff(response.flip) / (1/60));
        
        % how many interstimulus frames differed from the median?
        disp(sum(frames ~= median(frames)))
    end
    %
    
    
end

resp_trig_time = timing.trigger.flip_time(:);
timing.flip.channel = resp_trig_time;

tmp = timing.stimulus.flip_time(:);
tmp2 = timing.trigger.flip_time_02(:);

timing.flip.channel_02 = tmp(tmp2 > 0);




%% Compare to triggers from trigger channel
% Issues: it appears that the triggers were 'combined' with a lot
% crap in the trigger channels, making it impossible to work out the
% identity of all the triggers. More annoying is the fact that the timing
% of the triggers is not reliable as well....
% It seems that the trigger flip times are ahead of the triggers in the
% channels, to the extend these can be trusted. I.e. the trigger flips may
% be delayed compared to the trigger channel ??

timing.trigger.channel = all_trigger(meg_file, trig_chan);

%% Fiddle around with the photodiode
% Process trigger data:
% rescale to [0 1]

% Load the data:
phdata = sqdread(meg_file,'Channels',diode_chan);

phdata = phdata - min(phdata(:));
phdata = phdata / max(phdata(:));

% check whether triggers are indicated by a low value or a high value
if round(mean(phdata)) == 0, trigger_is_high = true;
else                     trigger_is_high = false; end

% if triggers are indicated by a low value, then invert the signal
if ~trigger_is_high, phdata = 1 - phdata; end

% threshold to binarize from analog recording
phdata = phdata > 0.75;


% differentiate to isolate trigger onsets, and pad to preserve length
phdata_d = [diff(phdata);0];

% rectify
phdata_d(phdata_d<0) = 0;

% mark the samples as 1 when any trigger channel signaled
any_trigger      = sum(phdata_d,2) > 0;

% list the sample numbers (time ponts) when any trigger channel signalled
any_trigger_inds = find(any_trigger);


timing.diode.channel = any_trigger_inds;


%% Find the triggers and photo diode flashes that correspond to eachother:
ph_trig_diff = bsxfun(@minus,timing.diode.channel,timing.trigger.channel(:,1)');

% ph_trig_diff = bsxfun(@minus, diode_flip_times2, any_trigger_inds');

% Select the flips and triggers that are closest to each other:
[min_diffs, min_idx] = min(abs(ph_trig_diff),[],1); % photo diode to trigger
[min_diffs2, min_idx2] = min(abs(ph_trig_diff),[],2); % Trigger to photo diode
% 
% idx = sub2ind(size(ph_trig_diff),min_idx,1:size(ph_trig_diff,2)); % used for obtaining the signed difference from ph_trig_diff
% idx2 = sub2ind(size(ph_trig_diff),1:size(ph_trig_diff,1),min_idx2');

% Remove any trigger that was more than 100 ms away from a photodiode
% flip, i.e. we simply did not send a photodiode flip for this
% trigger (i.e. they are from the init sequence):
keep_idx = min_idx(min_diffs < 100); % photo diode to trigger
keep_idx2 = min_idx2(min_diffs2 < 100); % trigger to photo diode
% 
% idx = idx(min_diffs < 100);
% idx2 = idx2(min_diffs2 < 100);
% min_diffs = min_diffs(min_diffs < 100);
% min_diffs2 = min_diffs2(min_diffs2 < 100);

keep_idx= unique(keep_idx);
% idx = idx(tmp);
% min_diffs = min_diffs(tmp);

keep_idx2 = unique(keep_idx2);
% idx2 = idx2(tmp);
% min_diffs2 = min_diffs2(tmp);
% These triggers were matched by the corresponding photodiode triggers:
timing.diode.diode2trigger = timing.diode.channel(keep_idx); % PHs with matched triggers

timing.trigger.trigger2diode = timing.trigger.channel(keep_idx2,:); % Triggers with matched PHs
timing.trigger.trigger2diode = timing.trigger.trigger2diode(diff([0;timing.trigger.trigger2diode(:,1)]) > 10,:);

timing.diode.diode2trigger_delay = timing.diode.diode2trigger - timing.trigger.trigger2diode(:,1);
timing.trigger.trigger2diode_delay = -timing.diode.diode2trigger_delay;


%% Find the triggers and flip times that correspond to eachother:
flip_trig_diff = bsxfun(@minus,resp_trig_time+timing.trigger.channel(1),timing.trigger.channel(:,1)');

% ph_trig_diff = bsxfun(@minus, diode_flip_times2, any_trigger_inds');

% Select the flips and triggers that are closest to each other:
[min_diffs, min_idx] = min(abs(flip_trig_diff),[],1); % Index to timing.diode.channel, with length of amount triggers, i.e. closest ph for every trigger
[min_diffs2, min_idx2] = min(abs(flip_trig_diff),[],2); % Index to trigger with length of amount PHs, i.e. closest trigger for every ph

% idx = sub2ind(size(flip_trig_diff),min_idx,1:size(flip_trig_diff,2)); % used for obtaining the signed difference from ph_trig_diff
% idx2 = sub2ind(size(flip_trig_diff),1:size(flip_trig_diff,1),min_idx2');

% Remove any trigger that was more than 100 ms away from a photodiode
% flip, i.e. we simply did not send a photodiode flip for this
% trigger (i.e. they are from the init sequence):
keep_idx = min_idx(min_diffs < 200); % PHs that are closer than 100 ms to the triggers for which they are the nearest
keep_idx2 = min_idx2(min_diffs2 < 200); % Triggers that are closer than 100 ms to the PHs for which thay are the nearest PH

keep_idx = unique(keep_idx);
keep_idx2 = unique(keep_idx2);
% These triggers were matched by the corresponding photodiode triggers:
timing.flip.flip2trigger = resp_trig_time(keep_idx)+timing.trigger.channel(1); % Flips with matched triggers
timing.trigger.trigger2flip = timing.trigger.channel(keep_idx2,:); % Triggers with matched flip times
timing.trigger.trigger2flip = timing.trigger.trigger2flip(diff([0;timing.trigger.trigger2flip(:,1)]) > 10,:);

timing.trigger.trigger2flip_delay = timing.flip.flip2trigger - timing.trigger.trigger2flip(:,1);
timing.flip.flip2trigger_delay = -timing.trigger.trigger2flip_delay;

%% Do the same with the triggers flip times and photo diode:
flip_ph_diff = bsxfun(@minus,resp_trig_time+timing.trigger.channel(1),timing.diode.channel');

% ph_trig_diff = bsxfun(@minus, diode_flip_times2, any_trigger_inds');

% Select the flips and triggers that are closest to each other:
[min_diffs, min_idx] = min(abs(flip_ph_diff),[],1); % Index to timing.diode.channel, with length of amount triggers, i.e. closest ph for every trigger
[min_diffs2, min_idx2] = min(abs(flip_ph_diff),[],2); % Index to trigger with length of amount PHs, i.e. closest trigger for every ph

% idx = sub2ind(size(flip_ph_diff),min_idx,1:size(flip_ph_diff,2)); % used for obtaining the signed difference from ph_trig_diff
% idx2 = sub2ind(size(flip_ph_diff),1:size(flip_ph_diff,1),min_idx2');

% Remove any trigger that was more than 100 ms away from a photodiode
% flip, i.e. we simply did not send a photodiode flip for this
% trigger (i.e. they are from the init sequence):
keep_idx = min_idx(min_diffs < 200); % PHs that are closer than 100 ms to the triggers for which they are the nearest
keep_idx2 = min_idx2(min_diffs2 < 200); % Triggers that are closer than 100 ms to the PHs for which thay are the nearest PH

keep_idx = unique(keep_idx);
keep_idx2 = unique(keep_idx2);
% These triggers were matched by the corresponding photodiode triggers:
timing.diode.diode2flip = timing.diode.channel(keep_idx2); % PHs with matched triggers

timing.flip.flip2diode = resp_trig_time(keep_idx)+timing.trigger.channel(1); % Triggers with matched flip times
timing.flip.flip2diode  = timing.flip.flip2diode (diff([0;timing.flip.flip2diode(:,1)]) > 10,:);


timing.flip.flip2diode_delay = timing.flip.flip2diode - timing.diode.diode2flip;
timing.diode.diode2flip_delay = -timing.flip.flip2diode;






%% Check stimulus timing/dropped frames etc
timing.stimulus.seqtime = timing.stimulus.seqtime(:);
timing.stimulus.flip_time = timing.stimulus.flip_time(:);

delay_sequence = timing.stimulus.flip_time - timing.stimulus.seqtime; % Difference between sequence and flips
d_delay_sequence = [0; diff(delay_sequence)]; % Differential, the flips that were delayed

if do_plot
    fh1 = figure;
    plot(delay_sequence)
    ylabel('Delay (milliseconds)');
    
    fh2 = figure;
    plot(timing.stimulus.flip_time(:), d_delay_sequence)
    xlabel('Stimulus frame flip time');
    ylabel('Added delay stimulus flips (ms)')
    
end
% 
% trigger_interval = diff(round(triggers(:,1)) ./ 1000);
% flip_interval = diff(timing.trigger.flip_time(:) ./ 1000);
% seq_interval = diff(timing.trigger.seq_time(:) ./ 1000);
% 
% % Amount of frames the triggers are lagging behind the flips:
% trigger.trigger_flip_delay = round((flip_interval - trigger_interval) ./ (1/params.display.frameRate));
% 
% % Amount of frames the flips are lagging behind the stimulus sequence
% trigger.flip_seq_delay = round((flip_interval - seq_interval) ./ (1/params.display.frameRate));
% 
% % Amount of frames the triggers are lagging behind the stimulus sequence
% trigger.trigger_seq_delay = round((trigger_interval - seq_interval) ./ (1/params.display.frameRate));

accum_delay = 0;
delay_idx = nan(1,10^5);
delay_idx_start = 1;
for n = 1:length(timing.stimulus.flip_time)
    
    cur_delay = timing.stimulus.flip_time(n) - timing.stimulus.seqtime(n);
    added_delay = round(cur_delay - accum_delay);
    accum_delay = accum_delay + added_delay;
    
    if added_delay > floor(1/params.display.frameRate * 1000);
        delay_start = round(timing.stimulus.seqtime(n) + (timing.stimulus.seqtime(n) - timing.stimulus.seqtime(n-1)))+1;
        delay_end =  delay_start + added_delay-1;
        
        delay_idx_end = delay_idx_start + added_delay-1;
        delay_idx(delay_idx_start:delay_idx_end) = delay_start : delay_end;
        
        
        delay_idx_start = delay_idx_start + added_delay;
        
    end
    
end

delay_idx = delay_idx(~isnan(delay_idx)); % indices of the time points that 
%were delayed. I.e. where the flip timing and sequence timing mismatch

if do_plot
    set(0,'CurrentFigure',fh2);
    hold on;
    plot(delay_idx, ones(size(delay_idx)) .* 100,'rx')
    
end

timing.delay = delay_sequence;
timing.d_delay = d_delay_sequence;
timing.delay_idx = delay_idx;





end























































