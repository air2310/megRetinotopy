function [PH_opt,VE_opt] = mprf_mostreliablephase(ft_data,opts,meg_resp)
% mprf_mostreliablephase - Determining the most reliable phase for
% every channel purely from the MEG data and the predicted response as the phase that gives highest variance explained for that channel
%
% input -
%       ft_data : Fourier transformed MEG data (4D data - Eg: 551 (frequency)x 140 (epochs) x 19 (repeats) x 157 (channels))
%       meg_resp : predicted MEG response (Eg: 140 (epochs) x 157 (channels))
%       opts : structure with information about, number of channels,
%              epochs, repeats, sampling rate, stimulus frequency, index of
%              stimulus frequency, metric
%

ph_range = -3.14:0.314:3.14; % range of values to search for the reference phase
VE_fit_alang = nan(size(ph_range,2),opts.n_chan);
ang_opt = nan(1,opts.n_chan);
PH_opt = nan(1,opts.n_chan);
VE_opt = nan(1,opts.n_chan);

for idx_ch = 1:opts.n_chan
    cur_ch = idx_ch;
    count_ang = 1;
    clear opt_ve VE_fit_opt;
    
    for idx_ang = ph_range
        tic;
        ref_ph = idx_ang;
        
        % Determine the phase reference amplitude for every channel
        [tseries_av, ~ , ~ ] = mprf_computemetric(ft_data(:,:,:,cur_ch),opts,ref_ph);
        
        % Compute the signed prediction
        % Fit the MEG predictions on the time series extracted above
        % Sizes (again...)
        n_bars = size(meg_resp{1},1);
        n_chan = size(meg_resp{1},2);
        n_roi = size(meg_resp{1},3);
        n_metric = size(tseries_av,3);
        
        fit_data_opt = nan(n_bars,n_chan);
        preds_opt =  nan(n_bars,n_chan);
        preds = nan(n_bars, 1, n_roi, n_metric);
        mean_ve = nan(1, n_roi, n_metric);
        fit_data = nan(n_bars,1,n_roi,n_metric);
        
        for this_chan = cur_ch
            for this_roi = 1:n_roi
                for this_metric = 1:n_metric
                    % Get MEG predictions:
                    cur_pred = meg_resp{1}(:,this_chan);
                    % Get current data:
                    cur_data = tseries_av(:,length(this_chan),this_metric);
                    % Exclude NANs
                    not_nan = ~isnan(cur_pred(:)) & ~isnan(cur_data(:));
                    
                    % Make X matrix (predictor)
                    if strcmpi(opts.metric, 'amplitude')
                        X = [ones(size(cur_pred(not_nan))) abs(cur_pred(not_nan))];
                    elseif strcmpi(opts.metric,'phase ref amplitude')
                        X = [ones(size(cur_pred(not_nan))) (cur_pred(not_nan))];
                    end
                    % Compute Beta's:
                    B = X \ cur_data(not_nan);
                    % Store the predicted times series:
                    preds(not_nan, idx_ch, this_roi, this_metric) =  X * B;
                    % Compute coefficient of determination (i.e. R square /
                    % variance explained):
                    cod_01 = 1- (var(cur_data(not_nan) - (X * B)) ./ var(cur_data(not_nan)));
                    mean_ve( idx_ch, this_roi, this_metric) = cod_01;
                    fit_data(:, idx_ch, this_roi, this_metric) = cur_data;
                end
            end
        end
        
        
        if notDefined('opt_ve')
            opt_ve = mean_ve(idx_ch);
        end
        if mean_ve(idx_ch) >= opt_ve
            ref_ph_opt = idx_ang;
            VE_fit_opt = mean_ve(idx_ch);
            count_opt = count_ang;
            opt_ve = VE_fit_opt;
            fit_data_opt(:,idx_ch) = fit_data(:, idx_ch, this_roi, this_metric);
            preds_opt(not_nan,idx_ch) = preds(not_nan,idx_ch);
        end
        if isnan(opt_ve)
            ref_ph_opt = NaN;
            VE_fit_opt = NaN;
            count_opt = NaN;
            opt_ve = NaN;
            fit_data_opt(:,idx_ch) = nan(n_bars,1);
            preds_opt(:,idx_ch) = nan(n_bars,1);
        end
        
        if count_ang == 1 && idx_ch == 1
            tmp_time = toc;
            tot_time = tmp_time * length(ph_range) * opts.n_looreps * n_chan;
            t_hms = datevec(tot_time./(60*60*24));
            fprintf('total time for the fitting : [%d %d %d %d %d %d] in Y M D H M S',t_hms);
            fprintf('\n');
        end
        VE_fit_alang(count_ang,idx_ch) = mean_ve(idx_ch);
        count_ang = count_ang+1;
    end
    
    ang_opt(idx_ch) = count_opt;
    PH_opt(idx_ch) = ref_ph_opt;
    VE_opt(idx_ch) = VE_fit_opt;
    
    % for plotting the phase values on polar cordinates
    %hold on;
    %x = (VE_fit_opt./max(VE_fit_opt)).*cos(ph_range);
    %y = (VE_fit_opt./max(VE_fit_opt)).*sin(ph_range);
    %plot(x,y,'ko-','MarkerSize',1);
    %hold on; plot([0 x(ang_opt(idx_ch))],[0 y(ang_opt(idx_ch))],'ro-');
    
    %heading = strcat('polarplot ch',num2str(idx_ch),{' '},num2str(rad2deg(PH_opt(idx_ch))));
    %set(figPoint_2_2,'Name',heading{1});
    %hold on;
    %x_bl = tseries_av_rel_amp_ch(periods.blank).* cos(tseries_av_rel_ph_ch(periods.blank));
    %y_bl = tseries_av_rel_amp_ch(periods.blank).* sin(tseries_av_rel_ph_ch(periods.blank));
    %plot(x_bl,y_bl,'ko');
    
    %hold off;
    %saveas(figPoint_2_2,strcat(save_dir,'/',num2str(idx_ch),'polarplot_ch'),'fig');
    %saveas(figPoint_2_2,strcat(save_dir,'/',num2str(idx_ch),'polarplot_ch',num2str(idx_ch)),'jpg');
    %saveas(figPoint_2_2,strcat(save_dir,'/',num2str(idx_ch),'polarplot_ch',num2str(idx_ch)),'eps');
    
end
end
