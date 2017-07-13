


stim_freq = 10;
harm_freq = stim_freq*2:stim_freq:stim_freq*10;
grid_freq = 60;
samp_rate = 1000;

sz = size(epoched_data.data);
n_time = sz(1);
n_bars = sz(2);
n_reps = sz(3);
n_chan = sz(4);


sl_idx = round(mprfFreq2Index(n_time,stim_freq,samp_rate));

stim_amp_av = nan(n_bars, n_chan);
stim_amp_std = nan(n_bars, n_chan);
stim_amp_ste = nan(n_bars, n_chan);

stim_co_av = nan(n_bars, n_chan);
stim_co_std = nan(n_bars, n_chan);
stim_co_ste = nan(n_bars, n_chan);

stim_ph_av = nan(n_bars, n_chan);
stim_ph_std = nan(n_bars, n_chan);
stim_ph_ste = nan(n_bars, n_chan);

fprintf('Processing %d stimulus periods:\n',sz(2));

for n = 1:n_bars
    fprintf('%d.',n)
    
    
    for nn = 1:n_chan
        cur_data = squeeze(epoched_data.data(:,n,:,nn));
        
         if any(~isnan(cur_data(:)))

             ft_data = fft(cur_data);
             
             ft_data = ft_data(1:1+fix(size(ft_data,1)/2),:);
             sc_amp = abs(ft_data);
             
             tmp_amp = 2*(sc_amp(sl_idx,:))/size(cur_data,1);
             
             stim_amp_av(n,nn) = nanmean(tmp_amp);
             stim_amp_std(n,nn) = nanstd(tmp_amp);
             stim_amp_ste(n,nn) = stim_amp_std(n,nn) ./ sqrt(sum(~isnan(tmp_amp)));
             
             sqrtsummagsq = sqrt(sum(sc_amp.^2));
             
             cur_sc_amp = sc_amp(sl_idx,:) ./ sqrtsummagsq;
             
             stim_co_av(n,nn) =  nanmean(cur_sc_amp);
             stim_co_std(n,nn) = nanstd(cur_sc_amp);
             stim_co_ste(n,nn) = stim_co_std(n,nn) ./ sqrt(sum(~isnan(cur_sc_amp)));
             
             
             cur_angle = angle(ft_data(sl_idx,:));
             
             stim_ph_av(n,nn) = nanmean(cur_angle);
             stim_ph_std(n,nn) = nanstd(cur_angle);
             stim_ph_ste(n,nn) = stim_ph_std(n,nn) ./ sqrt(sum(~isnan(cur_angle)));
             
             
             
             
         end
        
        
        
    end
    
    
    
end
fprintf('\n')


amp_r = nan(1,size(pred.meg_resp,2));

for n = 1:size(pred.meg_resp,2)
    
    sel = ~isnan(stim_amp_av(:,n)) & ~isnan(pred.meg_resp(:,n));
    mean(sel)
    
    tmp = corrcoef(zscore(stim_amp_av(sel,n)),zscore(abs(pred.meg_resp(sel,n))));
    amp_r(n) = tmp(1,2);
    
end

[amp_rs, s_idx] = sort(amp_r);

fh = figure;
megPlotMap(amp_r,[-1 1],fh, jet(256),'Correlation Prediction - Gain');

































