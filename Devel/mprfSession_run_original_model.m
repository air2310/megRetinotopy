function out = mprfSession_run_original_model(pred)

out = [];
if isnan(pred.model.params.n_iterations)
    pred.model.params.n_iterations = 1;
end

fh_sl = [];
fh_bb = [];
fh_sl_map = [];
fh_bb_map = [];

if ~exist('pred','var') || isempty(pred)
    pred = mprf__load_model_predictions;
    
end

model = pred.model;
meg_resp = pred.meg_resp;

opts.stim_freq = model.params.stim_freq;
opts.samp_rate = model.params.samp_rate;

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

sz = size(data.data);
opts.n_time = sz(1);
opts.n_bars = sz(2);
opts.n_reps = sz(3);
opts.n_chan = sz(4);


tmp = sum([model.params.do_sl model.params.do_bb]);
opts.metric = model.params.metric;
opts.idx = cell(1,tmp);

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


tseries_av = nan(opts.n_bars, opts.n_chan,size(opts.idx,2));
tseries_std = nan(opts.n_bars, opts.n_chan,size(opts.idx,2));
tseries_ste = nan(opts.n_bars, opts.n_chan,size(opts.idx,2));

if model.params.n_iterations > 1
    tseries_raw = nan(opts.n_bars, opts.n_reps, opts.n_chan,size(opts.idx,2));
end
fprintf('Processing %d stimulus periods:\n',opts.n_bars);

for this_metric = 1:size(opts.idx,2)
    cur_idx = opts.idx{this_metric};
    
    for this_bar = 1:opts.n_bars
        if mod(this_bar,10) == 0
            fprintf('%d.',this_bar)
        end
        for this_chan = 1:opts.n_chan
            cur_data = squeeze(ft_data(:,this_bar,:,this_chan));
            
            if strcmpi(opts.metric,'amplitude')
                
                if length(cur_idx) == 1
                    tmp =  squeeze(2*(abs(cur_data(cur_idx,:)))/opts.n_time);
                    n_nan = sum(isnan(tmp));
                    
                elseif length(cur_idx) > 1
                    tmp = 2*(abs(cur_data(cur_idx,:)))/opts.n_time;
                    tmp = squeeze(exp(nanmean(log(tmp.^2))));
                    n_nan = sum(isnan(tmp));
                    
                end
                
                if model.params.n_iterations > 1
                    tseries_raw(this_bar,:,this_chan,this_metric) = tmp;
                end
                
                tseries_av(this_bar,this_chan,this_metric) = nanmean(tmp);
                tseries_std(this_bar,this_chan,this_metric) = nanstd(tmp);
                tseries_ste(this_bar,this_chan,this_metric) = tseries_std(this_bar,this_chan,this_metric) ./ sqrt(n_nan);
                
            elseif strcmpi(opts.metric, 'coherence')
                error('Not implemented')
                
            elseif strcmpi(opts.metric, 'phase')
                error('Not implemented')
                
                
            else
                error('Not implemented')
            end
            
        end
    end
    fprintf('\n')
    
end

n_it = model.params.n_iterations;
n_cores = model.params.n_cores;

if  n_cores > 1
    if isempty(gcp('nocreate'))
        fprintf('No open pool found\n')
    else
        answer = questdlg('An open matlab pool is found. Do you want to close it or run on a single core',...
            'Open Matlab pool found','Close','Single core','Cancel','Close');
        
        switch lower(answer)
            
            case 'close'
                delete(gcp);
                
            case 'single core'
                pred.model.params.n_cores = 1;
                mprfSession_run_original_model(pred);
                
            case 'cancel'
                return
        end
    end
    
    
end


n_chan = size(meg_resp{1},2);
n_roi = size(meg_resp{1},3);
n_metric = size(tseries_av,3);


all_corr = nan(n_it, n_chan, n_roi, n_metric);

if n_cores > 1
    
    mpool = parpool(n_cores);
    
    
    parfor this_it = 1:n_it
        cur_idx = ceil(rand(1,opts.n_reps) .* opts.n_reps);
        
        
        for this_chan = 1:n_chan
            for this_roi = 1:n_roi
                for this_metric = 1:n_metric
                    if model.params.n_iterations > 1
                        cur_data = nanmean(tseries_raw(:,cur_idx,this_chan, this_metric),2);
                    
                    else
                        cur_data = tseries_av(:,this_chan,this_metric);
                    end
                    
                    cur_pred = meg_resp{1}(:,this_chan);
                    not_nan = ~isnan(cur_pred(:)) & ~isnan(cur_data(:));
                    
                    tmp = corrcoef(abs(cur_pred(not_nan)), cur_data(not_nan));
                    
                    all_corr(this_it, this_chan, this_roi, this_metric) = tmp(2);
                end
            end
        end
    end
    
    
    delete(mpool)
    
    
    

    
    
    
elseif n_cores == 1
    
    for this_it = 1:n_it
        for this_chan = 1:n_chan
            for this_roi = 1:n_roi
                for this_metric = 1:n_metric
                    
                    cur_pred = meg_resp{this_it}(:,this_chan);
                    cur_data = tseries_av(:,this_chan,this_metric);
                    not_nan = ~isnan(cur_pred(:)) & ~isnan(cur_data(:));
                    
                    tmp = corrcoef(abs(cur_pred(not_nan)), cur_data(not_nan));
                    
                    all_corr(this_it, this_chan, this_roi, this_metric) = tmp(2);
                end
            end
        end
    end

    
    
    
end


results.corr_mat = all_corr;
corr_ci = prctile(all_corr, [2.5 50 97.5],1);


if model.params.do_sl && model.params.do_bb
    results.corr_ci_sl = corr_ci(:,:,:,1);
    
    fh_sl = figure;
    plot(squeeze(results.corr_ci_sl(2,:,1)),'r','LineWidth',2);
    hold on;
    plot(1:157,squeeze(results.corr_ci_sl([1 3],:,1)),'k--')
    title('Stimulus locked')
    ylabel('Correlation')
    xlabel('Channel')
    
    results.corr_ci_bb = corr_ci(:,:,:,2);
    
    
    fh_bb = figure;
    plot(squeeze( results.corr_ci_bb(2,:,1)),'r','LineWidth',2);
    hold on;
    plot(1:157,squeeze(results.corr_ci_bb([1 3],:,1)),'k--')
    title('Broad band')
    ylabel('Correlation')
    xlabel('Channel')
    type = 'both';
    
    fh_sl_map = figure;
    megPlotMap(results.corr_ci_sl(2,:),[0 1],fh_sl_map,'jet','Correlation stimulus locked');
    
    
    fh_bb_map = figure;
    megPlotMap(results.corr_ci_bb(2,:),[0 1],fh_bb_map,'jet','Correlation broad band');
    
    
    
elseif ~model.params.do_sl && model.params.do_bb
    
    results.corr_ci_bb = corr_ci(:,:,:,2);
    
    
    fh_bb = figure;
    plot(squeeze( results.corr_ci_bb(2,:,1)),'r','LineWidth',2);
    hold on;
    plot(1:157,squeeze(results.corr_ci_bb([1 3],:,1)),'k--')
    title('Broad band')
    ylabel('Correlation')
    xlabel('Channel')
    
    type = 'Broadband';
    
        
    fh_bb_map = figure;
    megPlotMap(results.corr_ci_bb(2,:),[0 1],fh_bb_map,'jet','Correlation broad band');
    
    
elseif model.params.do_sl && ~model.params.do_bb
    results.corr_ci_sl = corr_ci(:,:,:,1);
    
    fh_sl = figure;
    plot(squeeze(results.corr_ci_sl(2,:,1)),'r','LineWidth',2);
    hold on;
    plot(1:157,squeeze(results.corr_ci_sl([1 3],:,1)),'k--')
    title('Stimulus locked')
    ylabel('Correlation')
    xlabel('Channel')
    
    type = 'Stimulus_locked';
    
    
    fh_sl_map = figure;
    megPlotMap(results.corr_ci_sl(2,:),[0 1],fh_sl_map,'jet','Correlation stimulus locked');
    
    
else
    
    
end


if isfield(pred,'cur_time') && ~isempty(pred.cur_time)
    cur_time = pred.cur_time;
    
else
   cur_time = mprf__get_cur_time; 
    
end



res_dir = mprf__get_directory('model_results');
main_dir = mprf__get_directory('main_dir');

if strcmpi(model.type,'run original model')
    rel_dir = 'original_model';
    
elseif strcmpi(model.type,'fix prf size')

   rel_dir = ['fixed_prf_size_' num2str(model.params.sigma_fix)];
    
end
save_dir = fullfile(main_dir, res_dir, rel_dir, ['Run_' type '_' cur_time]);
mkdir(save_dir);


if ~isempty(fh_sl)
   hgsave(fh_sl,fullfile(save_dir,'Stimulus_locked'));
end

if ~isempty(fh_bb)
   hgsave(fh_bb,fullfile(save_dir,'Broad_band'));
end

if ~isempty(fh_sl_map)
   hgsave(fh_sl_map,fullfile(save_dir,'Stimulus_locked_map'));
end

if ~isempty(fh_bb_map)
   hgsave(fh_bb_map,fullfile(save_dir,'Broad_band_map'));
end


save(fullfile(save_dir, 'Results'),'results','model')




% 
% if min(model.params.sigma_range) < 0;
%     one_idx = find(model.params.sigma_range == 0);
% else
%     one_idx = find(model.params.sigma_range == 1);
% end
% 
% [max_corr, midx] = max(all_corr,[],1);
% 
% max_corr = squeeze(max_corr);
% midx = squeeze(midx);
% 
% corr_at_one = squeeze(all_corr(one_idx,:,:,:));
% 
% corr_diff = abs(max_corr - corr_at_one);
% to_plot = 1;
% fh = figure;
% megPlotMap(corr_diff(:,to_plot),[0 ceil(max(corr_diff(:,to_plot)) ./ 0.05) .* 0.05],fh,'jet')
% 
% 
% fh1 = figure;
% megPlotMap(model.params.sigma_range(midx(:,to_plot)),...
%     [min(model.params.sigma_range) max(model.params.sigma_range)],fh1,'jet')
% 
% 
% fh2 = figure;
% megPlotMap(corr_at_one(:,to_plot),[0 1],fh2,'jet')
% 


end













