

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
first_flip_times = nan(1,nrepeats);
last_flip_times = nan(1,nrepeats);

flip_idx = cell(1,nrepeats);
stim_trig_idx = cell(1,nrepeats);
trig_flip_idx = cell(1,nrepeats);
trig_seq_time = cell(1,nrepeats);
trig_flip_time = cell(1,nrepeats);

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
            load(fullfile(param_dir,par_files(n).name),'response');
            
            first_flip_times(cond_match) = response.flip(1);
            last_flip_times(cond_match) = response.flip(end);
            
            flip_idx{cond_match} = round((response.flip - first_flip_times(1)) * 1000)'; % The times at which the screen flipped
            stim_trig_idx{cond_match} = round((stimulus.seqtiming+first_flip_times(cond_match) - first_flip_times(1))*1000); % The time for which a flip was requested
            % The code presents all frames in the sequence and does not
            % drop any of them, or does not speed the stimulus presentation
            % up to keep up with the timing. Therefore, missing the VBL
            % deadline will cause the flips to be delayed relative to the
            % stimulus sequence timing.
            trig_flip_idx{cond_match} = zeros(size(stim_trig_idx{cond_match}));
            trig_flip_idx{cond_match}(1:length(response.trig)) = response.trig';
            trig_seq_time{cond_match} = round((stimulus.seqtiming(stimulus.trigSeq > 0)+first_flip_times(cond_match) - first_flip_times(1))*1000);
            
            trig_flip_time{cond_match} = round((response.flip(stimulus.trigSeq > 0) - first_flip_times(1))*1000);
            
        end
    end
    
    
end

if sum(isnan(file_idx)) == 0
else
    error('Could not find a matching file for every repeat');
    
end

stim_02 = cell2mat(stim_trig_idx);
flip_02 = cell2mat(flip_idx);
trig_02 = cell2mat(trig_flip_idx);
trig_time = cell2mat(trig_seq_time);
trig_time_flip = cell2mat(trig_flip_time);

stim_03 = stim_02(:);
flip_03 = flip_02(:);
trig_03 = trig_02(:);
trig_time= trig_time(:);
trig_time_flip = trig_time_flip(:);

% figure; plot(stim_03, flip_03,'.')
% axis square
% xlabel('Frame time stimulus sequence (ms)');
% ylabel('Frame time flipes (ms)');

% aa = zeros(1,max(flip_03)+1);
% bb = aa;
% aa(stim_03+1) = 2;
% bb(flip_03+1) = 1;
% 
% figure; hold on;
% plot(aa);
% plot(bb,'r');

delay_sequence = flip_03 - stim_03;
figure; plot(delay_sequence)
ylabel('Delay (ms)');
xlabel('N frames since first frame')

some_diff_mat = [diff(round(triggers(:,1)) / 1000) diff(trig_time_flip / 1000) diff(trig_time / 1000)];
some_diff = [diff(round(triggers(:,1)) / 1000) - diff(trig_time / 1000)];

aa = zeros(size(delay_sequence));
aa(find(trig_03)) = [some_diff;0].*1000;

figure; hold on;
plot([0; diff(delay_sequence)])
plot(aa,'r');
ylabel('Delay (ms)');
xlabel('N frames since first frame')


return
%addpath(genpath('/Users/bpklein/Documents/PhD/pRF_and_MEG/denoiseproject'));

data_file = fullfile(raw_dir, fname);

flicker_freqs   = [30 20 15 12 10 6 5 4 ]; % Hz, i.e. images per second
ff_triggers = flicker_freqs;
blink_trigger = 3;
blank_trigger = 2;

info = sqdread(data_file,'info');

data_chan = 1:157;


n_epochs_total = length(triggers(:,1)) ./ 3 * (12 + 2 + 3);

epoch_id = nan(n_epochs_total,1);
first_epoch = zeros(size(epoch_id));
epoched_data = nan(1000, n_epochs_total, length(data_chan));
se_idx = nan(n_epochs_total, 2);


start_point = 1;
n = 0;
cnt = 0;
while start_point + 999 < info.SamplesAvailable && n < length(triggers(:,2))
    n = n + 1;
    cur_trigger = triggers(n,2);
    fprintf('Skipping %d ms prior to trigger nr %d (%d)\n',triggers(n,1) - start_point,...
        n, cur_trigger);
    
    start_point = triggers(n,1);
    
    
    if any(cur_trigger == ff_triggers)
        n_epochs = 12;
        
    elseif cur_trigger == blink_trigger
        n_epochs = 2;
    elseif cur_trigger == blank_trigger
        n_epochs = 3;
        
    else
        warning('Unknow trigger: %d, skipping',cur_trigger)
        continue
        
    end
    
    if n == length(triggers(:,1))
        end_point = start_point + 1000 * n_epochs;
    else
        end_point = triggers(n+1,1);
    end
    
    
    n_epochs_collected = 0;
    first_epoch_of_cur_trigger = true;
    
    while (start_point) + 999 < end_point && n_epochs_collected <= n_epochs
        cnt = cnt + 1;
        se_idx(cnt,:) = [start_point start_point + 1000];
        epoch_id(cnt) = cur_trigger;
        
        if first_epoch_of_cur_trigger
            first_epoch(cnt) = 1;
            first_epoch_of_cur_trigger = false;
        end
        
        tmp = sqdread(dfile,'Channels', data_chan, 'Samples', [start_point start_point + 999]);
        epoched_data(:,cnt,:) = permute(tmp,[1 3 2]);
        
        n_epochs_collected = n_epochs_collected + 1;
        start_point = start_point + 1000;
    end
end


epoched_data = epoched_data(:,1:cnt,:);
epoch_id = epoch_id(1:cnt);

save(fullfile(main_dir, cur_ddir,'Epoched_data'),'epoched_data');
save(fullfile(main_dir, cur_ddir,'Epoch_id'),'epoch_id');
save(fullfile(main_dir, cur_ddir,'First_epoch'),'first_epoch');






















return
trig_frame_seq = zeros(1,max(flip_03)+1);
trig_frame_seq(trig_time+1) = trig_03(trig_03 > 0);

trig_frame_flip = zeros(1,max(flip_03)+1);
trig_frame_flip(trig_time_flip+1) = trig_03(trig_03 > 0);

figure;hold on;
plot(trig_frame_seq)
plot(trig_frame_flip,'r');


figure; hold on;
plot(stim_02,flip_02,'k.')
plot([min(stim_02(:)) max(stim_02(:))], [min(flip_02(:)) max(flip_02(:))],'r')
axis square
xlabel('Stimulus sequence time');
ylabel('Flip time');


diff([stim_trig_idx{1}(1:10) flip_idx{1}(1:10)])
diff([stim_trig_idx{5}(1:10) flip_idx{5}(1:10)])


stimulus.seqtiming(1:10)







