function epoch = mprfEpochMEGData(raw_data, triggers, epoch_length, do_save, data_dir)

if ~exist('do_save','var') || isempty(do_save)
    do_save = false;
end

if do_save
    if ~exist('data_dir','var') || isempty(data_dir)
        data_dir = uigetdir(pwd,'Please select directory to store epoched data');
    end
end


n_epochs_total = length(triggers);
n_periods = 140;
n_repeats = n_epochs_total/n_periods;
n_channels = size(raw_data,1);

if n_repeats ~= round(n_repeats);
    warning('Amount of repeats is not a whole integer')
    n_repeats = round(n_repeats);
end

if isnumeric(epoch_length)
    n_data_points_per_epoch = epoch_length(2);
    skip_n_first_data_points = epoch_length(1);
    max_epoch_length = n_data_points_per_epoch;
    
elseif ischar(epoch_length)
    if strcmpi(epoch_length,'full')
        skip_n_first_data_points = 0;
        n_data_points_per_epoch = -1;
        
        max_epoch_length = max(max(diff(reshape(triggers(:,1),n_periods,n_repeats))));
        
        
    end
    
    
end

fprintf('Preallocating output variable ''epoched_data''.\n Estimated %d bites of memory.\n',...
    prod([max_epoch_length n_periods n_repeats n_channels 8]))

epoched_data = nan(max_epoch_length, n_periods, n_repeats, n_channels);

fprintf('Done\n')
epoch_id = nan(n_periods,1);
start_end_idx = nan(n_periods, n_repeats,2);

fprintf('Epoching data:\n');

fb = floor(length(triggers)/10:length(triggers)/10:length(triggers));


cur_step = 1;
for n = 1:length(triggers)
    
    cur_per = mod(n,n_periods);
    cur_rep = ceil(n/n_periods);
    
    if cur_per == 0
        cur_per = n_periods;
    end
    
    
    if mod(n,fb(cur_step)) == 0
        n_perc = n / length(triggers) * 100;
        fprintf('%0.2f%%\n',n_perc)
        cur_step = cur_step + 1;
    end
    
    if n <= n_periods
        epoch_id(n) = triggers(n,2);
    end
    
    start_idx = triggers(n,1)+skip_n_first_data_points;
    
    
    if n_data_points_per_epoch == -1
        if cur_per == n_periods % At the end of a period, take 1300 ms
            end_idx = start_idx + 1300-1;
            
        else
            end_idx = triggers(n+1,1)-1;
            
        end
        
    else
        end_idx = start_idx + n_data_points_per_epoch-1;
        
    end
    
    cur_ep_length = end_idx - start_idx +1;
    
    start_end_idx(cur_per,cur_rep,:) = [start_idx end_idx];
    epoched_data(1:cur_ep_length,cur_per,cur_rep,:) = permute(raw_data(:,start_idx:end_idx),[2 3 4 1]);
    
    
    
end
fprintf('Done\n')

cur_time = datestr(now);

cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';

epoch.data = epoched_data;
epoch.start_end = start_end_idx;
epoch.idx = epoch_id;

if do_save
    save(fullfile(data_dir, ['epoched_data_' cur_time '.mat']),'epoch','-v7.3');
end


end
% 
% n_channels = size(data_chan_data,1);
% 
% old_stim_id = 99999;
% delay_idx = timing.delay_idx + triggers(1); % delay idx starts at one, align it with MEG triggers
% 
% bl_period = nan((sum((triggers(:,2)) == 10) .* 1400),n_channels);
% stim_period = nan((sum(triggers(:,2) ~= 10 & triggers(:,2) ~=20) .* 1400),n_channels);
% 
% bl_st_idx= 1;
% stim_st_idx = 1;
% 
% fprintf('Processing data:\n')
% for n = 1:length(triggers)-1 % Skip last period for now, as we do not know where it ends
%     
%     
%     
%     
%     cur_id = triggers(n,2);
%     cur_st = triggers(n,1);
%     cur_end = triggers(n+1,1)-1;
%     
%     n_delays = ismember(cur_st:cur_end, delay_idx);
%     
%     if any(n_delays)
%         fprintf('Period with id %0.2f has %d delays, skipping\n',cur_id,sum(n_delays))
%         continue
%     end
%     
%     if cur_id == 20
%         fprintf('Skipping blink period\n');
%         
%     elseif cur_id == 10
%         tmp = data_chan_data(:,cur_st:cur_end)';
%         bl_end_idx = bl_st_idx + size(tmp,1)-1;
%         
%         bl_period(bl_st_idx:bl_end_idx,:) = bsxfun(@minus, tmp, mean(tmp));
%         
%         bl_st_idx = bl_end_idx+1;
%     elseif ismember(cur_id, 1:8)
%         
%         
%         if old_stim_id == cur_id % not the first stimulus id in a sequence
%             tmp = data_chan_data(:,cur_st:cur_end)';
%             stim_end_idx = stim_st_idx + size(tmp,1)-1;
%             
%             stim_period(stim_st_idx:stim_end_idx,:) = bsxfun(@minus, tmp, mean(tmp));
%             
%             stim_st_idx = stim_end_idx+1;
%         else % Possibly the first stimulus id in a sequence, i.e. the ones we want to skip
%             
%             fprintf('Skipping first bar position in sweep %d at %d\n',cur_id,cur_st);
%             
%         end
%         
%         old_stim_id = cur_id;
%         
%     else
%         warning('Unrecognized period id: %0.2f',cur_id);
%         
%         
%     end
%     
% end
% 
% fprintf('Done.\n');
% stim_period = stim_period(~isnan(stim_period(:,1)),:);
% bl_period = bl_period(~isnan(bl_period(:,1)),:);
% 
% 
% 
% 
% 
% 


















