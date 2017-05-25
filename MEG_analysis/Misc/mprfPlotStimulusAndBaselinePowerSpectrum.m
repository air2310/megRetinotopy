function mprfPlotStimulusAndBaselinePowerSpectrum(data_chan_data, triggers, timing)

n_channels = size(data_chan_data,1);

old_stim_id = 99999;
delay_idx = timing.delay_idx + triggers(1); % delay idx starts at one, align it with MEG triggers

bl_period = nan((sum((triggers(:,2)) == 10) .* 1400),n_channels);
stim_period = nan((sum(triggers(:,2) ~= 10 & triggers(:,2) ~=20) .* 1400),n_channels);

bl_st_idx= 1;
stim_st_idx = 1;

fprintf('Processing data:\n')
for n = 1:length(triggers)-1 % Skip last period for now, as we do not know where it ends
    
    cur_id = triggers(n,2);
    cur_st = triggers(n,1);
    cur_end = triggers(n+1,1)-1;
    
    n_delays = ismember(cur_st:cur_end, delay_idx);
    
    if any(n_delays)
        fprintf('Period with id %0.2f has %d delays, skipping\n',cur_id,sum(n_delays))
        continue
    end
    
    if cur_id == 20
        fprintf('Skipping blink period\n');
        
    elseif cur_id == 10
        tmp = data_chan_data(:,cur_st:cur_end)';
        bl_end_idx = bl_st_idx + size(tmp,1)-1;
        
        bl_period(bl_st_idx:bl_end_idx,:) = bsxfun(@minus, tmp, mean(tmp));        
        
        bl_st_idx = bl_end_idx+1;
    elseif ismember(cur_id, 1:8)
        
        
        if old_stim_id == cur_id % not the first stimulus id in a sequence
            tmp = data_chan_data(:,cur_st:cur_end)';
            stim_end_idx = stim_st_idx + size(tmp,1)-1;
            
            stim_period(stim_st_idx:stim_end_idx,:) = bsxfun(@minus, tmp, mean(tmp)); 
            
            stim_st_idx = stim_end_idx+1;
        else % Possibly the first stimulus id in a sequence, i.e. the ones we want to skip

            fprintf('Skipping first bar position in sweep %d at %d\n',cur_id,cur_st);
            
        end
        
        old_stim_id = cur_id;
        
    else
        warning('Unrecognized period id: %0.2f',cur_id);
        
        
    end
    
end

fprintf('Done.\n');
stim_period = stim_period(~isnan(stim_period(:,1)),:);
bl_period = bl_period(~isnan(bl_period(:,1)),:);


freq2idx = @(n_samples,F,s_rate) n_samples./(s_rate./F) +1;

stim_idx = freq2idx(size(stim_period,1),[10 13 60],1000); 
bl_idx = freq2idx(size(bl_period,1),[10 13 60],1000); 

freq_labels = 10:10:100;

stim_fs_axis = freq2idx(size(stim_period,1),freq_labels,1000);
bl_fs_axis = freq2idx(size(bl_period,1),freq_labels,1000);

ft_stim = fft(stim_period);
ft_stim = ft_stim(1:1+fix(size(stim_period,1)/2),:);

ft_bl = fft(bl_period);
ft_bl = ft_bl(1:1+fix(size(bl_period,1)/2),:);

power_spec_stim = (2 * abs(ft_stim) ./ size(stim_period,1)) .^2;
power_spec_bl = (2 * abs(ft_bl) ./ size(bl_period,1)) .^2;


if length(power_spec_stim) > round((max(stim_idx) .* 1.25))
   plot_idx = 1:round(max(stim_idx) .* 1.25);
else
    plot_idx = 1:length(power_spec_stim);
end

figure;hold on;
plot(log(power_spec_stim(plot_idx)))
ax_y_lim = ylim';
plot([stim_idx;stim_idx],ax_y_lim(:,ones(1, length(stim_idx))),'r--','LineWidth',2)
ylim(ax_y_lim);
set(gca,'Children',flipud(get(gca,'Children')),...
    'xtick',stim_fs_axis,...
    'xticklabel',freq_labels)
title('Stimulus period')
ylabel('Log power')
xlabel('Frequency')



if length(power_spec_bl) > round((max(bl_idx) .* 1.25))
   plot_idx = 1:round(max(bl_idx) .* 1.25);
else
    plot_idx = 1:length(power_spec_bl);
end

figure;hold on;
plot(log(power_spec_bl(plot_idx)))
ax_y_lim = ylim';
plot([bl_idx;bl_idx],ax_y_lim(:,ones(1, length(bl_idx))),'r--','LineWidth',2)
ylim(ax_y_lim);
set(gca,'Children',flipud(get(gca,'Children')),...
    'xtick',bl_fs_axis,...
    'xticklabel',freq_labels)
title('Baseline period')
ylabel('Log power')
xlabel('Frequency')

end