function data = mprf__do_reliability_analysis(model,type, data)
% Keep fpath and fname for storage, i.e. keep track of input data set.

if ~exist('type','var') || isempty(type)
    type = 'stimulus_locked';
    fprintf('No type defined, defaulting to stimulus locked\n')
end

do_sl = false;
do_bb = false;

fh_half =[];
fh_map = [];
fh_scans =[];

if strcmpi(type,'stimulus_locked')
    plot_title = 'Stimulus locked';
    do_sl = true;
elseif strcmpi(type, 'broadband')
    plot_title = 'Broadband';
    do_bb = true;
    
end

% For loading the MEG response to do the model based reference angle fit
if strcmpi(model.params.phase_fit,'model_fit')
    if ~exist('pred','var') || isempty(pred)
        pred = mprf__load_model_predictions;
        meg_resp = pred.meg_resp;
    end    
end

n_it = model.params.n_iterations_rel;
n_cores = model.params.n_cores;
opts.samp_rate = model.params.samp_rate;
opts.stim_freq = model.params.stim_freq;
do_split_half = model.params.reliability_split_half;
do_scans = model.params.reliability_scans;

if n_cores > 1
    if isempty(gcp('nocreate'))
        fprintf('No open pool found\n')
    else
        answer = questdlg('An open matlab pool is found. Do you want to close it or run on a single core',...
            'Open Matlab pool found','Close','Single core','Cancel','Close');
        
        switch lower(answer)
            
            case 'close'
                delete(gcp);
                
            case 'single core'
                model.params.n_cores = 1;
                mprf__do_reliability_analysis(model, type);
                
            case 'cancel'
                return
        end
    end
end


cur_dir = pwd;


if ~exist('data','var') || isempty(data)
    preproc_dir =  mprf__get_directory('meg_preproc');
    cd(preproc_dir)
    [fname, fpath] = uigetfile('*.mat', 'Select raw data to model');
    
    
    if fname == 0
        fprintf('No file selected, quitting\n');
        return
        
    end
    
    fprintf('Loading raw data...\n')
    tmp = load(fullfile(fpath, fname));
    var_name = fieldnames(tmp);
    data = tmp.(var_name{1});
    clear tmp
    
    cd(cur_dir)
end

periods.blank = [3:5 30:32 57:59 84:86 111:113 138:140];
periods.blink = [1 2 28 29 55 56 82 83 109 110 136 137];
periods.stim = setdiff(1:140,[periods.blink periods.blank]);
%upr = [73:81 115:123];
%lwr = [60:69 127:135];
%periods.stim = setdiff([upr lwr],[periods.blink periods.blank]);

sz = size(data.data);
n_time = sz(1);
n_bars = sz(2);
n_reps = sz(3);
n_chan = sz(4);
opts.n_time = n_time;
opts.metric = model.params.metric;
opts.n_bars = sz(2);
opts.n_reps = sz(3);
opts.n_chan = sz(4);

base_sel = false(1,n_bars);
base_sel(periods.stim) = true;

[stim_idx, bb_idx] = mprf__get_freq_indices(true, true, opts);

if model.params.do_sl && model.params.do_bb
    [opts.idx{1}, opts.idx{2}] = mprf__get_freq_indices(true, true, opts);
    
elseif model.params.do_sl && ~model.params.do_bb
    opts.idx{1} = mprf__get_freq_indices(true, false, opts);
elseif model.params.do_bb && ~model.params.do_sl
    [~, opts.idx{1}] = mprf__get_freq_indices(false, true, opts);
else
    error('Unknown option')
end

ft_data = mprf__fft_on_meg_data(data.data);

all_amp = nan(n_bars, n_reps, n_chan);

if strcmpi(model.params.metric,'amplitude')
    fprintf('Computing raw amplitudes\n')
    for this_chan = 1:n_chan
        for this_bar = 1:n_bars
            if do_sl
                all_amp(this_bar,:,this_chan) =  2*(abs(ft_data(stim_idx,this_bar,:,this_chan)))/n_time;
                
            elseif do_bb
                tmp = 2*(abs(ft_data(bb_idx,this_bar,:,this_chan)))/n_time;
                all_amp(this_bar,:,this_chan) =  exp(nanmean(log(tmp.^2)));
                
            end
            
        end
    end
elseif strcmpi(model.params.metric,'phase ref amplitude')
    % for phase ref amplitude/ for computing the most reliable phase per
    % channel
    opts.phs_metric = model.params.phase_fit;
    if strcmpi(opts.phs_metric,'data_fit')
        [PH_opt,VE_opt] = mprf_mostreliablephase_data(ft_data,opts,model);
    elseif strcmpi(opts.phs_metric,'model_fit')
        [PH_opt,VE_opt] = mprf_mostreliablephase(ft_data,meg_resp,opts,model);
    end
    
    fprintf('Computing phase referenced amplitudes\n')
    for this_chan = 1:n_chan
        for this_bar = 1:n_bars
            if do_sl
                tmp_phase_1 = angle(ft_data(stim_idx,this_bar,:,this_chan)); % take the phase at stimulus frequency, for all repeats. Should be 1 X 19
                %tmp2 = angle(nansum(exp(tmp(:)*1i))); % averages the phases across repeats, single value
                %if sum(isnan(tmp(:))) == opts.n_reps
                %    tmp2 = NaN;
                %end
                n_nan = sum(~isnan(tmp_phase_1));
                mst_rel_ang = PH_opt(this_chan);
                diff_ang = tmp_phase_1 - mst_rel_ang;
                diff_amp = cos(diff_ang);
                
                tmp_amp_1 = 2*(abs(ft_data(stim_idx,this_bar,:,this_chan)))/n_time;
                %tmp_all_amp(this_bar,:,this_chan) = tmp_amp_1;
                all_amp(this_bar,:,this_chan) = tmp_amp_1 .* diff_amp;
                
            elseif do_bb
                tmp = 2*(abs(ft_data(bb_idx,this_bar,:,this_chan)))/n_time;
                all_amp(this_bar,:,this_chan) =  exp(nanmean(log(tmp.^2)));
                
            end
            
        end
    end
end

warning('off','all')


if n_cores == 1
    
    if do_split_half
        fprintf('Starting split half test\n')
        
        fprintf('Doing %d iterations:\n',n_it)
        cor_stuff = nan(n_chan, n_it*2);
        
        n_01 = ceil(n_reps/2);
        
        fb = n_it .* (0.1:0.1:1) .*2;
        nfb = 0;
        for n = 1:2:(n_it*2)%model.params.n_iterations_rel
            
            if n > fb(1)
                nfb = nfb+1;
                fprintf('.%d %%.',10.*nfb)
                if length(fb) >1
                    fb = fb(2:end);
                end
                
            end
            
            tmp = randperm(n_reps);
            cur_idx_01 = tmp(1:n_01);
            cur_idx_02 = tmp(n_01+1:end);
            
            av_amp_01 = squeeze(nanmean(all_amp(:,cur_idx_01,:),2));
            av_amp_02 = squeeze(nanmean(all_amp(:,cur_idx_02,:),2));
            
            for this_chan = 1:n_chan
                sel = base_sel & ~isnan(av_amp_01(:,this_chan))' & ...
                    ~isnan(av_amp_02(:,this_chan))';
                
                
%                 X_01 = [ones(size(av_amp_01(sel, this_chan))) av_amp_01(sel, this_chan)];
%                 X_02 = [ones(size(av_amp_02(sel, this_chan))) av_amp_02(sel, this_chan)];                
%                 B_02 = X_01 \ X_02(:,2);
%                 B_01 = X_02 \ X_01(:,2);
%                 cod_02 = 1- (var(X_01(:,2) - (X_02 * B_02 )) ./ var(X_01(:,2)));
%                 cod_01 = 1- (var(X_02(:,2) - (X_01 * B_01)) ./ var(X_02(:,2)));
                

                 X_01 = [av_amp_01(sel, this_chan) ones(size(av_amp_01(sel, this_chan)))];
                 X_02 = [av_amp_02(sel, this_chan) ones(size(av_amp_02(sel, this_chan)))]; 
                 B_02 = polyfit(X_01(:,1),X_02(:,1),1);
                 B_01 = polyfit(X_02(:,1),X_01(:,1),1);
                 cod_02 = 1- (var(X_01(:,1) - (X_02 * B_02' )) ./ var(X_01(:,1)));
                 cod_01 = 1- (var(X_02(:,1) - (X_01 * B_01')) ./ var(X_02(:,1)));
                
                

                
                %tmp_r = corrcoef(av_amp_01(sel,this_chan), ...
                 %   av_amp_02(sel,this_chan));
                
                cor_stuff(this_chan,[n n+1]) = [cod_01 cod_02];
                
            end
            
        end
        
        results.split_half.corr_mat = cor_stuff;
        results.raw.amplitude = all_amp;

        med_corr = nanmedian(cor_stuff,2);
        
        figure;
        plot(med_corr)
        ylabel('Split-half correlation')
        xlabel('Channel')
        title(plot_title);
        
        fh = figure;
        megPlotMap(med_corr,[0 1],fh,'jet',['Split-half correlations ('  plot_title ')']);
        
        
    end
    
    if do_scans
        
%         
%         scan_corr = cell(1,n_reps-1);
%         all_av = squeeze(nanmean(all_amp,2));
%         scan_med_corr = nan(n_chan, n_reps-1);
%         
%         for nn = 1:n_reps
%             perms = mprf__get_unique_permutations(n_reps,nn,1000);
%             fprintf('N repetitions %d\n',nn);
%             fprintf('\nDoing %d iterations\n',length(perms))
%             
%             fb = length(perms) .* (0.1:0.1:1);
%             nfb = 0;
%             
%             
%             scan_corr{nn} = nan(n_chan,size(perms,1));
%             
%             for nnn = 1:size(perms,1)
%                 if nnn > fb(1)
%                     nfb = nfb+1;
%                     fprintf('.%d %%.',10.*nfb)
%                     if length(fb) >1
%                         fb = fb(2:end);
%                     end
%                     
%                 end
%                 
%                 
%                 cur_av = squeeze(nanmean(all_amp(:,perms(nnn,:),:),2));
%                 
%                 for this_chan = 1:n_chan
%                     sel = base_sel & ~isnan(cur_av(:,this_chan))' & ...
%                         ~isnan(all_av(:,this_chan))';
%                     tmp_r = corrcoef(cur_av(sel,this_chan), ...
%                         all_av(sel,this_chan));
%                     
%                     scan_corr{nn}(this_chan,nnn) = tmp_r(2);
%                     
%                 end
%                 
%                 
%             end
%             fprintf('Done \n')
%             
%             
%             
%         end
%         
%         
%         for nn = 1:n_reps
%             for  this_chan = 1:n_chan
%                 scan_med_corr(this_chan,nn) = median(scan_corr{nn}(this_chan,:));
%                 
%             end
%         end
%         
%         scan_med_corr_scaled = bsxfun(@times, scan_med_corr, med_corr);
%         
%         
%         
%         figure;
%         plot(scan_med_corr_scaled');
%         ylabel('Split-half correlation')
%         xlabel('N scans included')
%         
        
        
        scan_corr_02 = cell(1,n_reps-1);
        scan_med_corr_02 = nan(n_chan, n_reps-1);
        
        
        
        for nn = 2:n_reps
            perms = mprf__get_unique_permutations(n_reps,nn,1000);
            fprintf('N repetitions %d\n',nn);
            fprintf('\nDoing %d iterations\n',length(perms))
            
            fb = length(perms) .* (0.1:0.1:1) .*2;
            nfb = 0;
            
            
            scan_corr_02{nn} = nan(n_chan,size(perms,1).*2);
            
            for nnn = 1:2:(size(perms,1)*2)
                if nnn > fb(1)
                    nfb = nfb+1;
                    fprintf('.%d %%.',10.*nfb)
                    if length(fb) >1
                        fb = fb(2:end);
                    end
                    
                end
                n_01 = ceil(nn/2);
                
                cur_idx_01 = perms(ceil(nnn/2),1:n_01);
                cur_idx_02 = perms(ceil(nnn/2),n_01+1:end);
                
                av_amp_01 = squeeze(nanmean(all_amp(:,cur_idx_01,:),2));
                av_amp_02 = squeeze(nanmean(all_amp(:,cur_idx_02,:),2));
                
                for this_chan = 1:n_chan
                    sel = base_sel & ~isnan(av_amp_01(:,this_chan))' & ...
                        ~isnan(av_amp_02(:,this_chan))';
                    
                    
                    
%                 X_01 = [ones(size(av_amp_01(sel, this_chan))) av_amp_01(sel, this_chan)];
%                 X_02 = [ones(size(av_amp_02(sel, this_chan))) av_amp_02(sel, this_chan)];
%                 
%                 
%                 B_02 = X_01 \ X_02(:,2);
%                 B_01 = X_02 \ X_01(:,2);
%                 
%                 cod_02 = 1- (var(X_01(:,2) - (X_02 * B_02 )) ./ var(X_01(:,2)));
%                 cod_01 = 1- (var(X_02(:,2) - (X_01 * B_01)) ./ var(X_02(:,2)));
                
                 X_01 = [av_amp_01(sel, this_chan) ones(size(av_amp_01(sel, this_chan)))];
                 X_02 = [av_amp_02(sel, this_chan) ones(size(av_amp_02(sel, this_chan)))]; 
                 B_02 = polyfit(X_01(:,1),X_02(:,1),1);
                 B_01 = polyfit(X_02(:,1),X_01(:,1),1);
                 cod_02 = 1- (var(X_01(:,1) - (X_02 * B_02' )) ./ var(X_01(:,1)));
                 cod_01 = 1- (var(X_02(:,1) - (X_01 * B_01')) ./ var(X_02(:,1)));
                
                
                
                scan_corr_02{nn}(this_chan,[nnn nnn+1]) = [cod_01 cod_02];
                
                
                end
                
                
            end
            fprintf('Done \n')
            
            
            
        end
        
        results.scans.corr_mat = scan_corr_02;
        results.raw.amplitude = all_amp;

        for nn = 2:n_reps
            for  this_chan = 1:n_chan
                scan_med_corr_02(this_chan,nn-1) = median(scan_corr_02{nn}(this_chan,:));
                
            end
        end
        
        figure;
        plot(scan_med_corr_02')
        ylabel('Split-half correlation')
        xlabel('N scans')
        
        title(plot_title)
        
    end
    
elseif n_cores > 1
    
    mpool = parpool(n_cores);
    pctRunOnAll warning('off','all')

    
    cor_stuff = nan(n_chan, n_it);
    cor_stuff2 = nan(size(cor_stuff));
    
    n_01 = ceil(n_reps/2);
    if do_split_half
        fprintf('Starting split half test\n')
        fprintf('Doing %d iterations:\n',n_it)
        
        parfor n = 1:n_it%model.params.n_iterations_rel
            
            tmp = randperm(n_reps);
            cur_idx_01 = tmp(1:n_01);
            cur_idx_02 = tmp(n_01+1:end);
            
            av_amp_01 = squeeze(nanmean(all_amp(:,cur_idx_01,:),2));
            av_amp_02 = squeeze(nanmean(all_amp(:,cur_idx_02,:),2));
            
            for this_chan = 1:n_chan
                sel = base_sel & ~isnan(av_amp_01(:,this_chan))' & ...
                    ~isnan(av_amp_02(:,this_chan))';
              
%                 X_01 = [ones(size(av_amp_01(sel, this_chan))) av_amp_01(sel, this_chan)];
%                 X_02 = [ones(size(av_amp_02(sel, this_chan))) av_amp_02(sel, this_chan)];
%                 
%                 
%                 B_02 = X_01 \ X_02(:,2);
%                 B_01 = X_02 \ X_01(:,2);
%                 
%                 cod_02 = 1- (var(X_01(:,2) - (X_02 * B_02 )) ./ var(X_01(:,2)));
%                 cod_01 = 1- (var(X_02(:,2) - (X_01 * B_01)) ./ var(X_02(:,2)));

                 X_01 = [av_amp_01(sel, this_chan) ones(size(av_amp_01(sel, this_chan)))];
                 X_02 = [av_amp_02(sel, this_chan) ones(size(av_amp_02(sel, this_chan)))]; 
                 B_02 = polyfit(X_01(:,1),X_02(:,1),1);
                 B_01 = polyfit(X_02(:,1),X_01(:,1),1);
                 cod_02 = 1- (var(X_01(:,1) - (X_02 * B_02' )) ./ var(X_01(:,1)));
                 cod_01 = 1- (var(X_02(:,1) - (X_01 * B_01')) ./ var(X_02(:,1)));
                
                
                %tmp_r = corrcoef(av_amp_01(sel,this_chan), ...
                 %   av_amp_02(sel,this_chan));
                
                cor_stuff(this_chan,n) = cod_01;
                cor_stuff2(this_chan,n) = cod_02;
                
                
            end
            
        end
        cor_stuff = [cor_stuff cor_stuff2];
        results.split_half.corr_mat = cor_stuff;
        results.raw.amplitude = all_amp;
        med_corr = nanmedian(cor_stuff,2);
        
        fh_half = figure;
        plot(med_corr)
        ylabel('Split-half correlation')
        xlabel('Channel')
        title(plot_title)
        
        fh_map = figure;
        megPlotMap(med_corr,[0 1],fh_map,'jet', ['Split-half correlations ('  plot_title ')']);
        
       
    end
    
    if do_scans
%         fprintf('Starting scan reliability method 1\n')
%         scan_corr = cell(1,n_reps-1);
%         all_av = squeeze(nanmean(all_amp,2));
%         scan_med_corr = nan(n_chan, n_reps-1);
%         
%         
%         parfor nn = 1:n_reps
%             perms = mprf__get_unique_permutations(n_reps,nn,1000);
%             
%             fb = length(perms) .* (0.1:0.1:1);
%             nfb = 0;
%             
%             
%             scan_corr{nn} = nan(n_chan,size(perms,1));
%             
%             for nnn = 1:size(perms,1)
%                 if nnn > fb(1)
%                     nfb = nfb+1;
%                     fprintf('.%d %%.',10.*nfb)
%                     if length(fb) >1
%                         fb = fb(2:end);
%                     end
%                     
%                 end
%                 
%                 
%                 cur_av = squeeze(nanmean(all_amp(:,perms(nnn,:),:),2));
%                 
%                 for this_chan = 1:n_chan
%                     sel = base_sel & ~isnan(cur_av(:,this_chan))' & ...
%                         ~isnan(all_av(:,this_chan))';
%                     tmp_r = corrcoef(cur_av(sel,this_chan), ...
%                         all_av(sel,this_chan));
%                     
%                     scan_corr{nn}(this_chan,nnn) = tmp_r(2);
%                     
%                 end
%                 
%                 
%             end
%             fprintf('Done with %d iterations for %d repetitions\n', length(perms),nn)
%             
%         end
%         
%         
%         for nn = 1:n_reps
%             for  this_chan = 1:n_chan
%                 scan_med_corr(this_chan,nn) = median(scan_corr{nn}(this_chan,:));
%                 
%             end
%         end
%         
%         scan_med_corr_scaled = bsxfun(@times, scan_med_corr, med_corr);
%         
%         figure;
%         plot(scan_med_corr_scaled');
%         ylabel('Split-half correlation')
%         xlabel('N scans included')
%         
        
        
        scan_corr_02 = cell(1,n_reps-1);
        scan_med_corr_02 = nan(n_chan, n_reps-1);
                
        parfor nn = 2:n_reps
            perms = mprf__get_unique_permutations(n_reps,nn,1000);
            
            fb = length(perms) .* (0.1:0.1:1) .* 2;
            nfb = 0;
            
            
            scan_corr_02{nn} = nan(n_chan,size(perms,1).*2);
            
            for nnn = 1:2:(size(perms,1)*2)
                if nnn > fb(1)
                    nfb = nfb+1;
                    fprintf('.%d %%.',10.*nfb)
                    if length(fb) >1
                        fb = fb(2:end);
                    end
                    
                end
                n_01 = ceil(nn/2);
                
                cur_idx_01 = perms(ceil(nnn/2),1:n_01);
                cur_idx_02 = perms(ceil(nnn/2),n_01+1:end);
                
                av_amp_01 = squeeze(nanmean(all_amp(:,cur_idx_01,:),2));
                av_amp_02 = squeeze(nanmean(all_amp(:,cur_idx_02,:),2));
                
                for this_chan = 1:n_chan
                    sel = base_sel & ~isnan(av_amp_01(:,this_chan))' & ...
                        ~isnan(av_amp_02(:,this_chan))';
                    
%                     X_01 = [ones(size(av_amp_01(sel, this_chan))) av_amp_01(sel, this_chan)];
%                     X_02 = [ones(size(av_amp_02(sel, this_chan))) av_amp_02(sel, this_chan)];
%                     
%                     
%                     B_02 = X_01 \ X_02(:,2);
%                     B_01 = X_02 \ X_01(:,2);
%                     
%                     cod_02 = 1- (var(X_01(:,2) - (X_02 * B_02 )) ./ var(X_01(:,2)));
%                     cod_01 = 1- (var(X_02(:,2) - (X_01 * B_01)) ./ var(X_02(:,2)));
                    
                 X_01 = [av_amp_01(sel, this_chan) ones(size(av_amp_01(sel, this_chan)))];
                 X_02 = [av_amp_02(sel, this_chan) ones(size(av_amp_02(sel, this_chan)))]; 
                 B_02 = polyfit(X_01(:,1),X_02(:,1),1);
                 B_01 = polyfit(X_02(:,1),X_01(:,1),1);
                 cod_02 = 1- (var(X_01(:,1) - (X_02 * B_02' )) ./ var(X_01(:,1)));
                 cod_01 = 1- (var(X_02(:,1) - (X_01 * B_01')) ./ var(X_02(:,1)));
                  
                    
                    
                    scan_corr_02{nn}(this_chan,[nnn nnn+1]) = [cod_01 cod_02];
                    
                    
                    
                end
                
                
            end
            fprintf('Done with %d iterations for %d repetitions\n', length(perms),nn)
            
            
            
        end
        
        results.scans.corr_mat = scan_corr_02;
        results.raw.amplitude = all_amp;
        
        for nn = 2:n_reps
            for  this_chan = 1:n_chan
                scan_med_corr_02(this_chan,nn-1) = median(scan_corr_02{nn}(this_chan,:));
                
            end
        end
        
        fh_scans = figure;
        plot(scan_med_corr_02')
        ylabel('Split-half correlation')
        xlabel('N scans')
        title(plot_title)
    end
    
    pctRunOnAll warning('on','all')
    delete(mpool)
    
    
end
warning('on','all')

res_dir = mprf__get_directory('model_results');
main_dir = mprf__get_directory('main_dir');
rel_dir = 'reliability_checks';

cur_time = mprf__get_cur_time;
save_dir = fullfile(main_dir, res_dir, rel_dir, ['Run_' type '_' opts.phs_metric '_' cur_time]);
mkdir(save_dir);

if ~isempty(fh_half)
   hgsave(fh_half,fullfile(save_dir,'Split_half'));
end

if ~isempty(fh_map)
   hgsave(fh_map,fullfile(save_dir,'Split_half_scalp_map'));
end

if ~isempty(fh_scans)
   hgsave(fh_scans,fullfile(save_dir,'Split_half_scans'));
end

save(fullfile(save_dir, 'Results'),'results','model')


end