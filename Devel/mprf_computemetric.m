function [tseries_av, tseries_std, tseries_ste, tseries_av_rep] = mprf_computemetric(ft_data,opts,model,mst_rel_ang)
% Computes the metric
opts.n_bars = size(ft_data,2);
opts.n_reps = size(ft_data,3);
opts.n_chan = size(ft_data,4);

tseries_av = nan(opts.n_bars, opts.n_chan,size(opts.idx,2));
tseries_std = nan(opts.n_bars, opts.n_chan,size(opts.idx,2));
tseries_ste = nan(opts.n_bars, opts.n_chan,size(opts.idx,2));
tseries_av_rep = nan(opts.n_bars,opts.n_reps, opts.n_chan,size(opts.idx,2));
 for this_metric = 1:size(opts.idx,2)
        cur_idx = opts.idx{this_metric};
    % Loop over stimulus frames
    for this_bar = 1:opts.n_bars
        if mod(this_bar,10) == 0
            fprintf('%d.',this_bar)
        end
        % Loop over channels
        for this_chan = 1:opts.n_chan
%        for this_chan = mst_rel_channel    
            if ndims(ft_data)==4
               cur_data = squeeze(ft_data(:,this_bar,:,this_chan));
            elseif ndims(ft_data)==3
               cur_data = squeeze(ft_data(:,this_bar,:));
            end
            % If we want amplitude:
            if strcmpi(opts.metric,'amplitude')
                % If we want the amplitude from a single frequency
                % (i.e. stimulus locked)
                if length(cur_idx) == 1
                    tmp =  squeeze(2*(abs(cur_data(cur_idx,:)))/opts.n_time);
                    n_nan = sum(isnan(tmp));
                    % If we want the amplitude averaged across multiple
                    % frequencies (broad band):
                elseif length(cur_idx) > 1
                    tmp = 2*(abs(cur_data(cur_idx,:)))/opts.n_time;
                    tmp = squeeze(exp(nanmean(log(tmp.^2))));
                    n_nan = sum(isnan(tmp));
                    
                end
                % If we need the data for multiple iterations:
                if model.params.n_iterations > 1
                    tseries_raw(this_bar,:,this_chan,this_metric) = tmp;
                end
                
                % Store the average amplitude, it's standard deviation
                % and standard error across repetitions:
                tseries_av_rep(this_bar,:,this_chan,this_metric) = tmp;
                tseries_av(this_bar,this_chan,this_metric) = nanmean(tmp);
                tseries_std(this_bar,this_chan,this_metric) = nanstd(tmp);
                tseries_ste(this_bar,this_chan,this_metric) = tseries_std(this_bar,this_chan,this_metric) ./ sqrt(n_nan);
                
            elseif strcmpi(opts.metric, 'coherence')
                error('Not implemented')
                
            elseif strcmpi(opts.metric, 'phase')
                if length(cur_idx) == 1
                    tmp = angle(cur_data(cur_idx,:)); % take the phase at stimulus frequency, for all repeats. Should be 1 X 19
                    %tmp_1_2(this_bar,:,this_chan) = tmp;
                    tmp2 = angle(nansum(exp(tmp(:)*1i))); % averages the phases across repeats, single value
                    if sum(isnan(tmp(:))) == opts.n_reps
                        tmp2 = NaN;
                    end
                    n_nan = sum(~isnan(tmp));
                end
                
                % store the angular average and standard deviation of the phase
                tseries_av_rep(this_bar,:,this_chan,this_metric) = tmp;
                tseries_av(this_bar,this_chan,this_metric) = tmp2;
                tseries_std(this_bar,this_chan,this_metric) =  sqrt(-2*log((abs(nansum(exp(tmp(:).*1i))) ./ numel(n_nan))));
                % If you are implementing this, be sure to use the
                % circular/angular average and standard deviation:
                % ang_av = @(th) angle(sum(exp(th(:)*1i)));
                % ang_std = @(th) sqrt(-2*log((abs(sum(exp(th(:).*1i))) ./ numel(th))));
                % Look in mprf__coranal_on_meg_data.m for usage
                %
            elseif strcmpi(opts.metric,'phase ref amplitude')
                  % compute the phase of individual sensor and set it to
                  % the reference phase
                  if length(cur_idx) == 1
                    tmp = angle(cur_data(cur_idx,:)); % take the phase at stimulus frequency, for all repeats. Should be 1 X 19
                    tmp2 = angle(nansum(exp(tmp(:)*1i))); % averages the phases across repeats, single value
                    if sum(isnan(tmp(:))) == opts.n_reps
                        tmp2 = NaN;
                    end
                    n_nan = sum(~isnan(tmp));
                  end
                  diff_ang = tmp2 - mst_rel_ang; 
                  diff_amp = cos(diff_ang);
                  
                  % compute the new amplitude by considering the phase
                  if length(cur_idx) == 1
                      tmp_amp =  squeeze(2*(abs(cur_data(cur_idx,:)))/opts.n_time);
                      tmp2_amp = diff_amp .* tmp_amp;
                      n_nan = sum(isnan(tmp_amp));
                      % If we want the amplitude averaged across multiple
                      % frequencies (broad band):
                  elseif length(cur_idx) > 1
                      tmp_amp = 2*(abs(cur_data(cur_idx,:)))/opts.n_time;
                      tmp_amp = squeeze(exp(nanmean(log(tmp_amp.^2))));
                      tmp2_amp = diff_amp .* tmp_amp;
                      n_nan = sum(isnan(tmp_amp));
                      
                  end
                  % If we need the data for multiple iterations:
                  if model.params.n_iterations > 1
                      tseries_raw(this_bar,:,this_chan,this_metric) = tmp2_amp;
                  end
                  
                  % Store the average amplitude, it's standard deviation
                  % and standard error across repetitions:
                  tseries_av(this_bar,this_chan,this_metric) = nanmean(tmp2_amp);
                  tseries_std(this_bar,this_chan,this_metric) = nanstd(tmp2_amp);
                  tseries_ste(this_bar,this_chan,this_metric) = tseries_std(this_bar,this_chan,this_metric) ./ sqrt(n_nan);
                  
            else
                error('Not implemented')
            end
            
        end
    end
    fprintf('\n')
 end
end
