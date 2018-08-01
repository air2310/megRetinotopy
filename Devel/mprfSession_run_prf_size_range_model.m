function out = mprfSession_run_prf_size_range_model(pred)
out =[];

if ~exist('pred','var') || isempty(pred)
    pred = mprf__load_model_predictions;
    
end

if isnan(pred.model.params.n_iterations)
    pred.model.params.n_iterations = 1;
end

fh_sl = [];
fh_bb = [];
fh_sl_map = [];
fh_bb_map = [];


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


% for phase ref amplitude/ for computing the most reliable phase per
% channel
n_par_it = size(meg_resp,2);
n_it = model.params.n_iterations;
n_chan = size(meg_resp{1},2);
n_roi = size(meg_resp{1},3);
n_metric = size(tseries_av,3);

PH_opt = nan(n_par_it,n_chan);
VE_opt = nan(n_par_it,n_chan);
if strcmpi(opts.metric, 'phase ref amplitude')
    for this_par=1:n_par_it
        meg_resp_par{1} = meg_resp{this_par};
        [PH_opt(this_par,:),VE_opt(this_par,:)] = mprf_mostreliablephase(ft_data,meg_resp_par,opts,model);
    end
end


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
                
            elseif strcmpi(opts.metric, 'phase ref amplitude')
                
              for this_par =1:n_par_it
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
                mst_rel_ang = PH_opt(this_par,this_chan);
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
                    tseries_raw(this_bar,:,this_chan,this_metric,this_par) = tmp2_amp;
                end
                
                % Store the average amplitude, it's standard deviation
                % and standard error across repetitions:
                tseries_av(this_bar,this_chan,this_metric,this_par) = nanmean(tmp2_amp);
                tseries_std(this_bar,this_chan,this_metric,this_par) = nanstd(tmp2_amp);
                tseries_ste(this_bar,this_chan,this_metric,this_par) = tseries_std(this_bar,this_chan,this_metric,this_par) ./ sqrt(n_nan);
              end
            else
                error('Not implemented')
            end
            
        end
    end
    fprintf('\n')
    
end

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


all_corr = nan(n_it, n_par_it,n_chan, n_roi, n_metric);

if n_cores > 1
    
    mpool = parpool(n_cores);
    
    
    parfor this_it = 1:n_it
        cur_idx = ceil(rand(1,opts.n_reps) .* opts.n_reps);
        
        for this_par = 1:n_par_it
            for this_chan = 1:n_chan
                for this_roi = 1:n_roi
                    for this_metric = 1:n_metric
                        if model.params.n_iterations > 1
                            cur_data = nanmean(tseries_raw(:,cur_idx,this_chan, this_metric),2);
                            
                        else
                            cur_data = tseries_av(:,this_chan,this_metric);
                        end
                        
                        cur_pred = meg_resp{this_par}(:,this_chan);
                        not_nan = ~isnan(cur_pred(:)) & ~isnan(cur_data(:));
                        
                        tmp = corrcoef(abs(cur_pred(not_nan)), cur_data(not_nan));
                        
                        all_corr(this_it, this_par,this_chan, this_roi, this_metric) = tmp(2);
                    end
                end
            end
        end
    end
    
    delete(mpool)
    
    
    
    
    
    
    
elseif n_cores == 1
    
    for this_it = 1:n_it
        cur_idx = ceil(rand(1,opts.n_reps) .* opts.n_reps);
        
        
        for this_par = 1:n_par_it
            for this_chan = 1:n_chan
                for this_roi = 1:n_roi
                    for this_metric = 1:n_metric
                        
                        if model.params.n_iterations > 1
                            cur_data = nanmean(tseries_raw(:,cur_idx,this_chan, this_metric,this_par),2);
                            
                        else
                            cur_data = tseries_av(:,this_chan,this_metric,this_par);
                        end
                        
                        cur_pred = meg_resp{this_par}(:,this_chan);
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
                        preds(not_nan, this_chan, this_roi, this_metric) =  X * B;
                        % Compute coefficient of determination (i.e. R square /
                        % variance explained):
                        tmp = 1- (var(cur_data(not_nan) - (X * B)) ./ var(cur_data(not_nan)));
                        all_corr(this_it, this_par,this_chan, this_roi, this_metric) = tmp;
                        
                        
                        %tmp = corrcoef(abs(cur_pred(not_nan)), cur_data(not_nan));
                        %all_corr(this_it, this_par,this_chan, this_roi, this_metric) = tmp(2);
                    end
                end
            end
        end
    end
    
    
    
end

corr_ci = prctile(all_corr, [2.5 50 97.5],1);
[max_corr, mc_idx] = max(corr_ci(1,:,:,:,:),[],2);

if min(model.params.sigma_range) < 0;
    one_idx = find(model.params.sigma_range == 0);
else
    one_idx = find(model.params.sigma_range == 1);
end

corr_at_one = squeeze(corr_ci(2,one_idx,:,:,:));
range = [min(model.params.sigma_range) max(model.params.sigma_range)];


m_sigma_val = squeeze(model.params.sigma_range(mc_idx))';

results.corr_mat = all_corr;

if model.params.do_sl && model.params.do_bb
    results.corr_ci_sl = corr_ci(:,:,:,:,1);
    results.best_sigma_sl = m_sigma_val(:,1);
    results.best_corr_sl = squeeze(max_corr(:,:,:,:,1));
    results.corr_at_one_sl = corr_at_one(:,1);
    
    fh_sl = figure;
    plot(results.corr_at_one_sl,'r','LineWidth',2);
    hold on;
    plot(results.best_corr_sl,'k--')
    title('Stimulus locked')
    ylabel('Correlation')
    xlabel('Channel')
    
    
    fh_sl_map = figure;
    megPlotMap(results.best_sigma_sl,range,fh_sl_map,'jet','Best sigma difference stimulus locked',[],[],'interpmethod','nearest');

    
    results.corr_ci_bb = corr_ci(:,:,:,:,2);
    results.best_sigma_bb = m_sigma_val(:,2);
    results.best_corr_bb = squeeze(max_corr(:,:,:,:,2));
    results.corr_at_one_bb = corr_at_one(:,2);
    
    fh_bb = figure;
    plot(results.corr_at_one_bb,'r','LineWidth',2);
    hold on;
    plot(results.best_corr_bb,'k--')
    title('Stimulus locked')
    ylabel('Correlation')
    xlabel('Channel')
    
    type = 'both';
    
     
    
    fh_bb_map = figure;
    megPlotMap(results.best_sigma_bb,range,fh_bb_map,'jet','Best sigma difference broad band');
    
    
    
elseif ~model.params.do_sl && model.params.do_bb
    
    
    results.corr_ci_bb = corr_ci(:,:,:,:,2);
    results.best_sigma_bb = m_sigma_val(:,2);
    results.best_corr_bb = squeeze(max_corr(:,:,:,:,2));
    results.corr_at_one_bb = corr_at_one(:,2);
    
    fh_bb = figure;
    plot(results.corr_at_one_bb,'r','LineWidth',2);
    hold on;
    plot(results.best_corr_bb,'k--')
    title('Stimulus locked')
    ylabel('Correlation')
    xlabel('Channel')
    
    type = 'Broad_band';
    
     
    
    fh_bb_map = figure;
    megPlotMap(results.best_sigma_bb,range,fh_bb_map,'jet','Best sigma difference broad band');
    
    
elseif model.params.do_sl && ~model.params.do_bb
    results.corr_ci_sl = corr_ci(:,:,:,:,1);
    results.best_sigma_sl = m_sigma_val(:,1);
    results.best_corr_sl = squeeze(max_corr(:,:,:,:,1));
    results.corr_at_one_sl = corr_at_one(:,1);
    
    fh_sl = figure;
    plot(results.corr_at_one_sl,'r','LineWidth',2);
    hold on;
    plot(results.best_corr_sl,'k--')
    title('Stimulus locked')
    ylabel('Correlation')
    xlabel('Channel')
    
    type = 'Stimulus_locked';
    
    fh_sl_map = figure;
    megPlotMap(results.best_sigma_sl,range,fh_sl_map,'jet','Best sigma difference stimulus locked');

    
%     corr_tmp=squeeze(all_corr);
%     for tmp_cnt=1:n_par_it
%         fh_sl_map=figure;
%         megPlotMap(corr_tmp(tmp_cnt,:),[0 0.6],fh_sl_map,'jet',...
%             'Phase ref fit stimulus locked',[],[],'interpmethod','nearest');
%         
%         saveas(fh_sl_map,strcat('sl_map_',num2str(tmp_cnt),'.jpg'));
%     end
    
    load(which('meg160_example_hdr.mat'))
    layout = ft_prepare_layout([],hdr);
    xpos = layout.pos(1:157,1);
    ypos = layout.pos(1:157,2);
    %chan_occ = find(ypos<1 & xpos<1);
    chan_occ = find(ypos<0 & xpos<1);
    % fh_allchan = mprfPlotHeadLayout(chan_occ);
    % saveas(fh_allchan,'Occchan_hp.jpg');
    
    tmp_corr = nan(n_par_it,size(chan_occ,1));
    for i=1:n_par_it
        tmp_corr(i,:) = squeeze(all_corr(1,i,chan_occ));
    end
    corr_avg = nanmean(tmp_corr,2);
    corr_std = nanstd(tmp_corr,0,2);
    corr_ste = corr_std ./ sqrt(size(tmp_corr,2));
    corr_CI = 1.96 .* corr_ste;
    corr_ci_occ = prctile(tmp_corr,[5 50 95],2);
    
    fh_sl_occ = figure;
    %plot(1:n_par_it,corr_ci_occ(:,2),'r','Linewidth',3);
    plot(1:n_par_it,corr_avg,'r','Linewidth',3);
    %hold on; stem(1:n_par_it,[corr_ci_occ(:,1),corr_ci_occ(:,3)],'b')
    %hold on; stem(1:n_par_it,corr_avg,'k');
    %lb = corr_ci_occ(:,2) - corr_ci_occ(:,1);
    %ub = -corr_ci_occ(:,2) + corr_ci_occ(:,3);
    %hold on; errorbar(1:n_par_it,corr_ci_occ(:,2),lb,ub);
    hold on; errorbar(1:n_par_it,corr_avg,corr_CI);
    grid on;
    ylim([0.05 0.14]);
    set(gca,'XTick',1:n_par_it);
    set(gca,'XTickLabel',model.params.sigma_range);
    title('Variance explained per sigma ratio');
    xlabel('Sigma ratio');
    ylabel('Variance explained');
    
else
    
    
end


if isfield(pred,'cur_time') && ~isempty(pred.cur_time)
    cur_time = pred.cur_time;
    
else
    cur_time = mprf__get_cur_time;
    
end



res_dir = mprf__get_directory('model_results');
main_dir = mprf__get_directory('main_dir');
rel_dir = 'prf_size_range';
 
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

if ~isempty(fh_sl_occ)
    hgsave(fh_sl_occ,fullfile(save_dir,'Stimulus_locked_VEvsSig_occ'));
end


save(fullfile(save_dir, 'Results'),'results','model')

end
