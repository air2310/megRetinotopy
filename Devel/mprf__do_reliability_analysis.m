
% Keep fpath and fname for storage, i.e. keep track of input data set.
global mprfSESSION

preproc_dir =  mprf__get_directory('meg_preproc');
cur_dir = pwd;

cd(preproc_dir)
[fname, fpath] = uigetfile('*.mat', 'Select raw data to model');
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

n_01 = ceil(n_reps/2);

opts.do_av = true;
opts.do_std = false;
opts.do_ste = false;
opts.do_ft = false;

if strcmpi(model.params.metric, 'amplitude')
    tmp_meas_02 = '_amp';
    
else
    error('Metric not implemented')
    
end

if model.params.do_sl
    tmp_meas_01 = 'stim';
    
end

opts.measure = [tmp_meas_01 tmp_meas_02];

samp_rate = 1000;
stim_freq = 10;
opts.stim_idx = round(mprfFreq2Index(n_time,stim_freq,samp_rate));

n_it = 1000;

av_amp_01 = nan(n_bars, n_chan, n_it);
av_amp_02 = nan(n_bars, n_chan, n_it);


base_sel = false(1,n_bars);
base_sel(periods.stim) = true;

all_amp = nan(n_bars, n_reps, n_chan);

for this_chnl = 1:n_chan
    for this_bar = 1:n_bars
        tmp = fft(data.data(:,this_bar,:,this_chnl));
        tmp = tmp(1:1+fix(n_time/2),:,:,:);
        all_amp(this_bar,:,this_chnl) =  2*(abs(tmp(opts.stim_idx,:,:)))/n_time;
        
        
    end
end

fprintf('Doing %d iterations:\n',n_it)
cor_stuff = nan(n_chan, n_it);
    

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



return

if model.params.reliability_split_half
    for n = 1:n_it%model.params.n_iterations_rel
        fprintf('.%d.',n)
        for this_chnl = 1:n_chan
            for this_bar = 1:n_bars
                tmp = randperm(n_reps);
                cur_idx_01 = tmp(1:n_01);
                cur_idx_02 = tmp(n_01+1:end);
                
                data_01.ft = squeeze(ft_data(:,this_bar,cur_idx_01,this_chnl));
                data_02.ft = squeeze(ft_data(:,this_bar,cur_idx_02,this_chnl));
                data_01.sc_amp = abs(data_01.ft);
                data_02.sc_amp = abs(data_02.ft);
                data_01.size_01 = size_01;
                data_02.size_01 = size_01;
                
                it_01 = mprf__compute_meg_measure(data_01,opts);
                it_02 = mprf__compute_meg_measure(data_02,opts);
                
                av_amp_01(this_bar,this_chnl,n) = it_01.av;
                av_amp_02(this_bar,this_chnl,n) = it_02.av;
                
                
                
            end
            
            sel = base_sel & ~isnan(av_amp_01(:,this_chnl,n))' & ...
                ~isnan(av_amp_02(:,this_chnl,n))';
            
            tmp_r = corrcoef(av_amp_01(sel,this_chnl,n), ...
                av_amp_02(sel,this_chnl,n));
            cor_stuff(this_chnl,n) = tmp_r(2);
        end
        
        
    end
end

fprintf('\n')

ch_to_plot = 15;
it_to_plot = 1;




if model.params.reliability_scans
    
    
    
    
    
end


















