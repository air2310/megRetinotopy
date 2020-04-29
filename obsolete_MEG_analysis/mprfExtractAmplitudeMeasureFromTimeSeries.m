function data_out = mprfExtractAmplitudeMeasureFromTimeSeries(data, measure)


fprintf('Computing measure: %s\n',measure)
if strcmpi(measure, 'amp_ratio')
    
    %% Preprocess data:
    
stim_freq = 10;
harm_freq = stim_freq*2:stim_freq:stim_freq*10;
grid_freq = 60;

bl_freq = [9 11];
bl_freq_02 = [9 10 11];
exl_freq = [harm_freq grid_freq];
all_idx = 1:1+fix(size(data.data,1)/2);

include_idx = setdiff(all_idx, round(mprfFreq2Index(size(data.data,1),exl_freq,1000)));

nb_idx_tmp = round(mprfFreq2Index(size(data.data,1),bl_freq,1000));
sl_idx_tmp = round(mprfFreq2Index(size(data.data,1),stim_freq,1000));
%bl_idx_tmp = round(mprfFreq2Index(size(data.data,1),bl_freq_02,1000));



sz = size(data.data);
snr_tmp = nan(sz([2 4]));
snr_median = nan(sz([2 4]));
snr_data = nan(sz([2 4 3]));
sl_co = nan(sz([2 4]));
sl_amp = nan(sz([2 4]));
sl_co_std = nan(sz([2 4]));
sl_co_ste = nan(sz([2 4]));
sl_ph = nan(sz([2 4]));
sl_amp_std = nan(sz([2 4]));
sl_amp_ste = nan(sz([2 4]));
sl_ph_std = nan(sz([2 4]));
sl_ph_ste = nan(sz([2 4]));

fprintf('Processing %d stimulus periods:\n',sz(2));
for n = 1:sz(2)
    fprintf('%d.',n)
    
    for nn = 1:sz(4);
        cur_data = squeeze(data.data(:,n,: , nn));
        
        if any(~isnan(cur_data(:)))
            data_av_tmp = nanmean(cur_data,2);
            ft_av_tmp = fft(data_av_tmp);
            ft_data = fft(cur_data);
            
            amp_av_tmp = (2 .* abs(ft_av_tmp) ./ size(ft_av_tmp,1));
            amp_data = (2 .* abs(ft_data) ./ size(ft_data,1));
            
            snr_tmp(n,nn) = exp(nanmean(log(amp_av_tmp(sl_idx_tmp)))) ...
                ./ exp(nanmean(log(amp_av_tmp(nb_idx_tmp))));
            
            snr_data(n,nn,:) = exp(nanmean(log(amp_data(sl_idx_tmp,:)))) ...
                ./ exp(nanmean(log(amp_data(nb_idx_tmp,:))));
            
            
            spec_amp_nb = abs(amp_data(nb_idx_tmp,:,:));
            spec_amp_sl = abs(amp_data(sl_idx_tmp,:,:));
            
            av_nb = squeeze(exp(nanmean(log(spec_amp_nb))));
            av_sl = squeeze(spec_amp_sl);
            
            snr_median(n,nn) = nanmedian(av_sl,2) ./ nanmedian(av_nb,2);
            
            
            %%% CorAnal:
            ft = fft(cur_data);
            ft = ft(1:1+fix(size(ft, 1)/2), :);
            
            % This quantity is proportional to the amplitude
            %
            scaledAmp = abs(ft);
            
            % This is in fact, the correct amplitude
            %
            cur_amp = 2*(scaledAmp(sl_idx_tmp,:))/size(cur_data,1);
            sl_amp(n,nn) = nanmean(cur_amp);
            sl_amp_std(n,nn) = nanstd(cur_amp);
            sl_amp_ste(n,nn) = sl_amp_std(n,nn) ./ sqrt(sum(~isnan(cur_amp)));
            
            sqrtsummagsq = sqrt(sum(scaledAmp(include_idx,:).^2));
            
            cur_sc_amp = scaledAmp(sl_idx_tmp,:) ./ sqrtsummagsq;
            
            sl_co(n,nn) =  nanmean(cur_sc_amp);
            sl_co_std(n,nn) = nanstd(cur_sc_amp);
            sl_co_ste(n,nn) = sl_co_std(n,nn) ./ sqrt(sum(~isnan(cur_sc_amp)));
             
            
            cur_angle = angle(ft(sl_idx_tmp,:));
            
            sl_ph(n,nn) = nanmean(cur_angle);
            sl_ph_std(n,nn) = nanstd(cur_angle);
            sl_ph_ste(n,nn) = sl_ph_std(n,nn) ./ sqrt(sum(~isnan(cur_angle)));
            
            
            
        end
    end
end
fprintf('\n')

data_out.snr_tmp = snr_tmp;
data_out.snr_data = snr_data;
data_out.snr_median = snr_median;

data_out.sl_co = sl_co;
data_out.sl_co_std = sl_co_std;
data_out.sl_co_ste = sl_co_ste;

data_out.sl_ph = sl_ph;
data_out.sl_ph_std = sl_ph_std;
data_out.sl_ph_ste = sl_ph_ste;

data_out.sl_amp = sl_amp;
data_out.sl_amp_ste = sl_amp_ste;
data_out.sl_amp_std = sl_amp_std;

else
    error('Unrecognized measure type')
    
    
    
    
end

fprintf('Done.\n')























