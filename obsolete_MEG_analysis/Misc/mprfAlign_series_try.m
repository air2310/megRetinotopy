%#ok<*BDSCI>

epoch_length = 1000;
epoch_ps_length = 100;
epoch_skip = 151;
epoch_idx = epoch_skip : epoch_skip + epoch_length-1;
n_channels = size(data_chan_data,1);

n_repeats = 19;
n_periods = size(triggers,1) ./ n_repeats;

epoched_ps = nan(n_repeats,n_periods,epoch_ps_length+1,n_channels);
delay_idx = timing.delay_idx + triggers(1); % delay idx starts at one, align it with MEG triggers

fprintf('Processing data:\n')
for n = 1:length(triggers)-1 % Skip last period for now, as we do not know where it ends
    
    
    cur_per = mod(n,n_periods);
    cur_rep = ceil(n/n_periods);
    
    if cur_per == 0
        cur_per = n_periods;
    end
    
    cur_id = triggers(n,2);
    cur_idx = triggers(n,1) + epoch_idx;
    
    %     cur_st = triggers(n,1);
    %     cur_end = triggers(n+1,1)-1;
    
    n_delays = ismember(cur_idx, delay_idx);
    
    if any(n_delays)
        fprintf('Period with id %0.2f has %d delays, skipping\n',cur_id,sum(n_delays))
        continue
    end
    
    tmp = data_chan_data(:,cur_idx)';
    tmp_ft = fft(tmp,[],1);
    epoched_ps(cur_rep,cur_per,:,:) = log((2*abs(tmp_ft(1:epoch_ps_length+1,:)) ./ size(tmp,1)).^2);
    
end

fprintf('Done.\n');

base = 51; 
end_pl = base+9;
for n = base:end_pl
    figure;
    plot(squeeze(epoched_ps(5,15,:,n)))
    title(sprintf('%d',n));
end


some_idx = [31 52 53 54 56 57 52 48 39 34 36 37 38 24 27 15 11 19 33];
figure;
plot(squeeze(epoched_ps(1,6,1:30,some_idx)))

return
t1 = 15;
t2 = 13;
cur_fs = 10;

freq2idx = @(n_samples,F,s_rate) n_samples./(s_rate./F) +1;


st_01 = triggers(t1,1)+151;
st_02 = triggers(t2,1)+151;

end_01 = triggers(t1,1)+1150;
end_02 = triggers(t2,1)+1150;

ep_01 = data_chan_data(1,st_01:end_01);
ep_02 = data_chan_data(1,st_02:end_02);


cur_idx_01 = round(freq2idx(size(ep_01,2),cur_fs,1000));
cur_idx_02 = round(freq2idx(size(ep_02,2),cur_fs,1000));


shift = mod(diff([st_01 st_02]),100);

shift(shift > 50) = shift(shift > 50) - 100;
ep_02_sh = mprfShiftTimeSeries(ep_02,1000,shift);

ft_01 = fft(ep_01);
ft_02 = fft(ep_02);
ft_02_sh = fft(ep_02_sh);

ft_01 = ft_01(1:floor(size(ep_01,2)/2)+1);
ft_02 = ft_02(1:floor(size(ep_02,2)/2)+1);
ft_02_sh = ft_02_sh(1:floor(size(ep_02_sh,2)/2)+1);

ph_01 = angle(ft_01(cur_idx_01));
ph_02 = angle(ft_02(cur_idx_02));
ph_02_sh = angle(ft_02_sh(cur_idx_02));


[ph_01 ph_02 ph_02_sh]


return


