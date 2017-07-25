function mprf__do_reliability_analysis(model)
% Keep fpath and fname for storage, i.e. keep track of input data set.

n_it = model.params.n_iterations_rel;
n_cores = model.params.n_cores;
samp_rate = model.params.samp_rate;
stim_freq = model.params.stim_freq;
do_split_half = model.params.reliability_split_half;
do_scans = model.params.reliability_scans;


preproc_dir =  mprf__get_directory('meg_preproc');
cur_dir = pwd;

cd(preproc_dir)
[fname, fpath] = uigetfile('*.mat', 'Select raw data to model');
fprintf('Loading raw data...\n')
tmp = load(fullfile(fpath, fname));
var_name = fieldnames(tmp);
data = tmp.(var_name{1});
clear tmp

cd(cur_dir)


periods.blank = [3:5 30:32 57:59 84:86 111:113 138:140];
periods.blink = [1 2 28 29 55 56 82 83 109 110 136 137];
periods.stim = setdiff(1:140,[periods.blink periods.blank]);


sz = size(data.data);
n_time = sz(1);
n_bars = sz(2);
n_reps = sz(3);
n_chan = sz(4);

stim_idx = round(mprfFreq2Index(n_time,stim_freq,samp_rate));


base_sel = false(1,n_bars);
base_sel(periods.stim) = true;

all_amp = nan(n_bars, n_reps, n_chan);
fprintf('Computing raw amplitudes\n')
for this_chnl = 1:n_chan
    for this_bar = 1:n_bars
        tmp = fft(data.data(:,this_bar,:,this_chnl));
        tmp = tmp(1:1+fix(n_time/2),:,:,:);
        all_amp(this_bar,:,this_chnl) =  2*(abs(tmp(stim_idx,:,:)))/n_time;
        
        
    end
end

if n_cores == 1
    
    if do_split_half
        fprintf('Starting split half test\n')
        
        fprintf('Doing %d iterations:\n',n_it)
        cor_stuff = nan(n_chan, n_it);
        
        n_01 = ceil(n_reps/2);
        
        fb = n_it .* (0.1:0.1:1);
        nfb = 0;
        for n = 1:n_it%model.params.n_iterations_rel
            
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
            
            for this_chnl = 1:n_chan
                sel = base_sel & ~isnan(av_amp_01(:,this_chnl))' & ...
                    ~isnan(av_amp_02(:,this_chnl))';
                tmp_r = corrcoef(av_amp_01(sel,this_chnl), ...
                    av_amp_02(sel,this_chnl));
                
                cor_stuff(this_chnl,n) = tmp_r(2);
                
            end
            
        end
        
        med_corr = nanmedian(cor_stuff,2);
        
        figure;
        plot(med_corr)
        ylabel('Split-half correlation')
        xlabel('Channel')
        
        fh = figure;
        megPlotMap(med_corr,[0 1],fh,'jet','Split-half correlations')
        
        
    end
    
    if do_scans
        
        
        scan_corr = cell(1,n_reps-1);
        all_av = squeeze(nanmean(all_amp,2));
        scan_med_corr = nan(n_chan, n_reps-1);
        
        for nn = 1:n_reps
            perms = mprf__get_unique_permutations(n_reps,nn,1000);
            fprintf('N repetitions %d\n',nn);
            fprintf('\nDoing %d iterations\n',length(perms))
            
            fb = length(perms) .* (0.1:0.1:1);
            nfb = 0;
            
            
            scan_corr{nn} = nan(n_chan,size(perms,1));
            
            for nnn = 1:size(perms,1)
                if nnn > fb(1)
                    nfb = nfb+1;
                    fprintf('.%d %%.',10.*nfb)
                    if length(fb) >1
                        fb = fb(2:end);
                    end
                    
                end
                
                
                cur_av = squeeze(nanmean(all_amp(:,perms(nnn,:),:),2));
                
                for this_chnl = 1:n_chan
                    sel = base_sel & ~isnan(cur_av(:,this_chnl))' & ...
                        ~isnan(all_av(:,this_chnl))';
                    tmp_r = corrcoef(cur_av(sel,this_chnl), ...
                        all_av(sel,this_chnl));
                    
                    scan_corr{nn}(this_chnl,nnn) = tmp_r(2);
                    
                end
                
                
            end
            fprintf('Done \n')
            
            
            
        end
        
        
        for nn = 1:n_reps
            for  this_chnl = 1:n_chan
                scan_med_corr(this_chnl,nn) = median(scan_corr{nn}(this_chnl,:));
                
            end
        end
        
        scan_med_corr_scaled = bsxfun(@times, scan_med_corr, med_corr);
        

        
        figure;
        plot(scan_med_corr_scaled');
        ylabel('Split-half correlation')
        xlabel('N scans included')
        
        
        
        scan_corr_02 = cell(1,n_reps-1);
        scan_med_corr_02 = nan(n_chan, n_reps-1);
        
        
        
        for nn = 2:n_reps
            perms = mprf__get_unique_permutations(n_reps,nn,1000);
            fprintf('N repetitions %d\n',nn);
            fprintf('\nDoing %d iterations\n',length(perms))
            
            fb = length(perms) .* (0.1:0.1:1);
            nfb = 0;
            
            
            scan_corr_02{nn} = nan(n_chan,size(perms,1));
            
            for nnn = 1:size(perms,1)
                if nnn > fb(1)
                    nfb = nfb+1;
                    fprintf('.%d %%.',10.*nfb)
                    if length(fb) >1
                        fb = fb(2:end);
                    end
                    
                end
                n_01 = ceil(nn/2);
                
                cur_idx_01 = perms(nnn,1:n_01);
                cur_idx_02 = perms(nnn,n_01+1:end);
                
                av_amp_01 = squeeze(nanmean(all_amp(:,cur_idx_01,:),2));
                av_amp_02 = squeeze(nanmean(all_amp(:,cur_idx_02,:),2));
                
                for this_chnl = 1:n_chan
                    sel = base_sel & ~isnan(av_amp_01(:,this_chnl))' & ...
                        ~isnan(av_amp_02(:,this_chnl))';
                    tmp_r = corrcoef(av_amp_01(sel,this_chnl), ...
                        av_amp_02(sel,this_chnl));
                    
                    
                    scan_corr_02{nn}(this_chnl,nnn) = tmp_r(2);
                    
                    
                end
                
                
            end
            fprintf('Done \n')
            
            
            
        end
        
        for nn = 2:n_reps
            for  this_chnl = 1:n_chan
                scan_med_corr_02(this_chnl,nn-1) = median(scan_corr_02{nn}(this_chnl,:));
                
            end
        end
        
        figure;
        plot(scan_med_corr_02')
        
    end
    
elseif n_cores > 1
    
    mpool = parpool(n_cores);
    
    
    cor_stuff = nan(n_chan, n_it);
    
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
            
            for this_chnl = 1:n_chan
                sel = base_sel & ~isnan(av_amp_01(:,this_chnl))' & ...
                    ~isnan(av_amp_02(:,this_chnl))';
                tmp_r = corrcoef(av_amp_01(sel,this_chnl), ...
                    av_amp_02(sel,this_chnl));
                
                cor_stuff(this_chnl,n) = tmp_r(2);
                
            end
            
        end
        
        med_corr = nanmedian(cor_stuff,2);
                
        figure;
        plot(med_corr)
        ylabel('Split-half correlation')
        xlabel('Channel')
        
        fh = figure;
        megPlotMap(med_corr,[0 1],fh,'jet','Split-half correlations')
        
        
    end
    
    if do_scans
        fprintf('Starting scan reliability method 1\n')
        scan_corr = cell(1,n_reps-1);
        all_av = squeeze(nanmean(all_amp,2));
        scan_med_corr = nan(n_chan, n_reps-1);
        
        
        parfor nn = 1:n_reps
            perms = mprf__get_unique_permutations(n_reps,nn,1000);
                
            fb = length(perms) .* (0.1:0.1:1);
            nfb = 0;
            
            
            scan_corr{nn} = nan(n_chan,size(perms,1));
            
            for nnn = 1:size(perms,1)
                if nnn > fb(1)
                    nfb = nfb+1;
                    fprintf('.%d %%.',10.*nfb)
                    if length(fb) >1
                        fb = fb(2:end);
                    end
                    
                end
                
                
                cur_av = squeeze(nanmean(all_amp(:,perms(nnn,:),:),2));
                
                for this_chnl = 1:n_chan
                    sel = base_sel & ~isnan(cur_av(:,this_chnl))' & ...
                        ~isnan(all_av(:,this_chnl))';
                    tmp_r = corrcoef(cur_av(sel,this_chnl), ...
                        all_av(sel,this_chnl));
                    
                    scan_corr{nn}(this_chnl,nnn) = tmp_r(2);
                    
                end
                
                
            end
            fprintf('Done with %d iterations for %d repetitions\n', length(perms),nn)

        end
        
        
        for nn = 1:n_reps
            for  this_chnl = 1:n_chan
                scan_med_corr(this_chnl,nn) = median(scan_corr{nn}(this_chnl,:));
                
            end
        end
        
        scan_med_corr_scaled = bsxfun(@times, scan_med_corr, med_corr);

        figure;
        plot(scan_med_corr_scaled');
        ylabel('Split-half correlation')
        xlabel('N scans included')
        
        
        
        scan_corr_02 = cell(1,n_reps-1);
        scan_med_corr_02 = nan(n_chan, n_reps-1);
        
        
        fprintf('Starting scan reliability method 2\n')
        parfor nn = 2:n_reps
            perms = mprf__get_unique_permutations(n_reps,nn,1000);
             
            fb = length(perms) .* (0.1:0.1:1);
            nfb = 0;
            
            
            scan_corr_02{nn} = nan(n_chan,size(perms,1));
            
            for nnn = 1:size(perms,1)
                if nnn > fb(1)
                    nfb = nfb+1;
                    fprintf('.%d %%.',10.*nfb)
                    if length(fb) >1
                        fb = fb(2:end);
                    end
                    
                end
                n_01 = ceil(nn/2);
                
                cur_idx_01 = perms(nnn,1:n_01);
                cur_idx_02 = perms(nnn,n_01+1:end);
                
                av_amp_01 = squeeze(nanmean(all_amp(:,cur_idx_01,:),2));
                av_amp_02 = squeeze(nanmean(all_amp(:,cur_idx_02,:),2));
                
                for this_chnl = 1:n_chan
                    sel = base_sel & ~isnan(av_amp_01(:,this_chnl))' & ...
                        ~isnan(av_amp_02(:,this_chnl))';
                    tmp_r = corrcoef(av_amp_01(sel,this_chnl), ...
                        av_amp_02(sel,this_chnl));
                    
                    
                    scan_corr_02{nn}(this_chnl,nnn) = tmp_r(2);
                    
                    
                end
                
                
            end
            fprintf('Done with %d iterations for %d repetitions\n', length(perms),nn)
            
            
            
        end
    end
    delete(mpool)
    
    
    for nn = 2:n_reps
        for  this_chnl = 1:n_chan
            scan_med_corr_02(this_chnl,nn-1) = median(scan_corr_02{nn}(this_chnl,:));
            
        end
    end
    
    figure;
    plot(scan_med_corr_02')
    
end




end