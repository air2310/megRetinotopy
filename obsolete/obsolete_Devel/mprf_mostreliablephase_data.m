function [PH_opt,VE_opt] = mprf_mostreliablephase_data(ft_data,opts)
% mprf_mostreliablephase_data - Determining the most reliable phase for
% every channel purely from the MEG data without using the predicted response
%
% input -
%       ft_data : Fourier transformed MEG data (4D data - Eg: 551 (frequency)x 140 (epochs) x 19 (repeats) x 157 (channels))
%       opts : structure with information about, number of channels,
%              epochs, repeats, sampling rate, stimulus frequency, index of
%              stimulus frequency, metric
%

ang_opt = nan(1,opts.n_chan);
PH_opt = nan(1,opts.n_chan);
VE_opt = nan(1,opts.n_chan);
ph_range = -3.14:0.0314:3.14;
VE_all = nan(length(ph_range),opts.n_chan);

for idx_ch = 1:opts.n_chan
    
    count_ang = 1;
    tseries_av_ph_ref_amp_var = nan(size(ph_range));
    clear opt_var VE_fit ref_ph_opt;
    for idx_ang = ph_range
        
        % Reference phase
        ref_ph = idx_ang;
        
        % Phase referenced amplitude
        [tseries_av_ph_ref_amp, ~,~] = mprf_computemetric(ft_data(:,:,:,idx_ch),opts,ref_ph);
        
        % variance of the phase referenced amplitude values for a given channel with a given reference phase
        tseries_av_ph_ref_amp_var(:,count_ang) = nanvar(tseries_av_ph_ref_amp);
        
        % Set the reference phase as the one with highest variance
        if notDefined('opt_var')
            opt_var = tseries_av_ph_ref_amp_var(:,count_ang);
        end
        if tseries_av_ph_ref_amp_var(:,count_ang) >= opt_var
            ref_ph_opt = idx_ang;
            VE_fit_opt = tseries_av_ph_ref_amp_var(:,count_ang);
            opt_var = tseries_av_ph_ref_amp_var(:,count_ang);
            count_opt = count_ang;
        end
        count_ang = count_ang+1;
    end
    VE_all(:,idx_ch) = tseries_av_ph_ref_amp_var'; % each point on the black ellipse
    ang_opt(idx_ch) = count_opt;
    PH_opt(idx_ch) = ref_ph_opt; % Reference phase for every channel
    VE_opt(idx_ch) = VE_fit_opt;
    
end



end
