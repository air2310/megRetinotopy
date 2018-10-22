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
fh_sl_VEparams = [];
fh_bb_VEparams = [];

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
opts.n_iter = model.params.n_iterations;

% params for the phase fitting
opts.phs_metric = model.params.phase_fit; % for choosing between model or data based fit
phase_fit_loo.do = model.params.phase_fit_loo;
phase_fit_loo.type = model.params.loo_type; %'amp' or 'VE';
phfittype = '';
opts.n_looreps = 1;

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

% Separate the ft_data to train and test data set to do a leave one out
% analysis
idx_train = nan(opts.n_reps,opts.n_reps-1);
idx_test = nan(opts.n_reps,1);
if phase_fit_loo.do == 1
    opts.n_looreps = opts.n_reps;
    for this_loorep = 1:opts.n_looreps
        idx_test(this_loorep,1) = this_loorep;
        idx_train(this_loorep,:) = setdiff(1:opts.n_reps,idx_test(this_loorep));
    end
else
    idx_test = 1:opts.n_reps;
    idx_train = 1:opts.n_reps;
end

% starting the parallel pool for running par for loop for reference
% phase fitting for Leave one out repetitions
n_cores = model.params.n_cores;
if n_cores > 1
    mpool = parpool(n_cores);
    pctRunOnAll warning('off','all');
end


if model.params.n_iterations > 1 || n_cores > 1
    tseries_raw = nan(opts.n_bars, opts.n_reps, opts.n_chan,size(opts.idx,2));
end

fprintf('Processing %d stimulus periods:\n',opts.n_bars);


% for phase ref amplitude/ for computing the most reliable phase per
% channel
n_par_it = size(meg_resp,2);
n_it = model.params.n_iterations;
n_chan = size(meg_resp{1},2);
n_roi = size(meg_resp{1},3);
n_metric = size(opts.idx,2);

%%
tic;
if strcmpi(opts.metric, 'phase ref amplitude')
    PH_opt_tmp = nan(opts.n_looreps,n_chan,n_par_it);
    VE_opt_tmp = nan(opts.n_looreps,n_chan,n_par_it);
    if strcmpi(opts.phs_metric,'data_fit')
        % Determining the reference phase from MEG data alone, as the
        % phase that gives maximum variance for phase referenced
        % amplitudes for 140 epochs for that particular channel.
        
        [PH_opt_tmp,VE_opt_tmp] = mprf_mostreliablephase_data(ft_data,opts);
        
        PH_opt = repmat(PH_opt_tmp,1,1,n_par_it);
        VE_opt = repmat(VE_opt_tmp,1,1,n_par_it);
    
    elseif strcmpi(opts.phs_metric,'model_fit')
        % Determing the reference phase from MEG data and the predicted
        % repsonses, as the phase that gives highest variance explained
        % for a particular channel.
        % A reference phase is computed for every pRF size range value and
        % for every leave one out repetition value.
        
        % Leave one out computation of reference phase
        % Step 1: Take 18 repeats (leaving the first repeat out). Calculate
        %         the reference phase for every channel (157) as the angle
        %         that gives highest variance explained (COD between the phase referenced amplitude metric
        %         and the predicted MEG response for every pRF size range value)
        % Step 2: Repeat this 19 times for every leave one out condition.
        % Step 3: Repeat this for all the pRF size range scaling values (here it is 19 ratios)
        
        for this_par=1:n_par_it % For every pRf size range value
            meg_resp_par{1} = meg_resp{this_par}; % There is a meg prediction array of 140x157 for every pRF size range value. 
                                                  % One is taken at a time to compute the reference phase.
            
            parfor this_loorep = 1:opts.n_looreps % For every leave one out condition
                
                [PH_opt_tmp(this_loorep,:,this_par),VE_opt_tmp(this_loorep,:,this_par)] = mprf_mostreliablephase(ft_data(:,:,idx_train(this_loorep,:),:),opts,meg_resp_par);
            end    
                if phase_fit_loo.do == 1
                    phfittype = 'lo';
%                     for cur_chan =1:opts.n_chan
%                         xbins = -3.14:0.314:3.14-0.314;
%                         [N,x] = hist(PH_opt_tmp(:,cur_chan,this_par),xbins);
%                         idx_max = find(N == max(N));
%                         idx_close = (x(idx_max(1))-0.1745 < PH_opt_tmp(:,cur_chan,this_par)') & (PH_opt_tmp(:,cur_chan)' < x(idx_max(1))+0.1745);
%                         PH_opt_tmp(~idx_close,cur_chan,this_par) = wrapToPi(PH_opt_tmp(~idx_close,cur_chan,this_par) + 3.14);
%                     end
                end
                
            %end
            
        end
        
    end
    results.PH_opt = PH_opt_tmp;
    results.VE_opt = VE_opt_tmp;
end
tot_time = toc;
t_hms = datevec(tot_time./(60*60*24));
fprintf('total time for the fitting : [%d %d %d %d %d %d] in Y M D H M S',t_hms);

% Ending the parallel pool
if isempty(gcp('nocreate'))
    fprintf('?? no open pool found ??');
else
    delete(gcp);
end
%%

% to check something for the leave one out condition, its better to load
% the precomputed reference phases.
%  load('modeling/results/prf_size_range/Run_Stimulus_locked_Model_fit_lo_25_Sep_2018_17_39_01/Results.mat');
%  PH_opt_tmp = results.PH_opt;

tseries_av = nan(opts.n_bars, opts.n_chan,size(opts.idx,2),n_par_it,opts.n_looreps);
tseries_std = nan(opts.n_bars, opts.n_chan,size(opts.idx,2),n_par_it,opts.n_looreps);
tseries_ste = nan(opts.n_bars, opts.n_chan,size(opts.idx,2),n_par_it,opts.n_looreps);
  
for this_loorep = 1:opts.n_looreps % For leave one out computation of the phases
    if strcmpi(opts.metric, 'phase ref amplitude') && strcmpi(opts.phs_metric,'model_fit')
        PH_opt = PH_opt_tmp(this_loorep,:,:);
    end
    
    for this_par = 1:n_par_it
        [tseries_av(:,:,:,this_par,this_loorep), tseries_std(:,:,:,this_par,this_loorep) ,...
            tseries_ste(:,:,:,this_par,this_loorep)] = mprf_computemetric(ft_data(:,:,idx_test...
            (this_loorep,:),:),opts,PH_opt(1,:,this_par));
    end
end
if strcmpi(phase_fit_loo.type,'Amp')
    tseries_av = nanmean(tseries_av,5);
    opts.n_looreps = size(tseries_av,5);
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
                            cur_data = nanmean(tseries_raw(:,cur_idx,this_chan, this_metric,this_par),2);
                            
                        else
                            if strcmpi(opts.metric,'phase ref amplitude')
                                cur_data = tseries_av(:,this_chan,this_metric,this_par);
                            else
                                cur_data = tseries_av(:,this_chan,this_metric);
                            end
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
                        %preds(not_nan, this_chan, this_roi, this_metric) =  X * B;
                        % Compute coefficient of determination (i.e. R square /
                        % variance explained):
                        tmp = 1- (var(cur_data(not_nan) - (X * B)) ./ var(cur_data(not_nan)));
                        all_corr(this_it, this_par,this_chan, this_roi, this_metric) = tmp;
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
                            if strcmpi(opts.metric,'phase ref amplitude')
                                cur_data = tseries_av(:,this_chan,this_metric,this_par);
                            else
                                cur_data = tseries_av(:,this_chan,this_metric);
                            end
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

range = [min(model.params.x0_range) max(model.params.x0_range)];

m_sigma_val = squeeze(model.params.sigma_range(mc_idx));

corr_at_one = squeeze(corr_ci(2,one_idx,:,:,:));


results.corr_mat = all_corr;

load(which('meg160_example_hdr.mat'))
layout = ft_prepare_layout([],hdr);
xpos = layout.pos(1:157,1);
ypos = layout.pos(1:157,2);
%chan_occ = find(ypos<1 & xpos<1);
chan_occ = find(ypos<0 & xpos<1);

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
    megPlotMap(results.best_sigma_sl,range,fh_sl_map,'jet','Best sigma difference stimulus locked');
    
    
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
    
    
    % Position ratio vs Variance explained for sl
    tmp_corr_sl = nan(n_par_it,size(chan_occ,1));
    for i=1:n_par_it
        tmp_corr_sl(i,:) = squeeze(all_corr(:,i,chan_occ,:,1));
    end
    corr_avg_sl = nanmean(tmp_corr_sl,2);
    corr_std_sl = nanstd(tmp_corr_sl,0,2);
    corr_ste_sl = corr_std_sl ./ sqrt(size(tmp_corr_sl,2));
    corr_CI_sl = 1.96 .* corr_ste_sl;
    fh_sl_VEparams = figure;
    plot(1:n_par_it,corr_avg_sl,'r','Linewidth',3);
    hold on; errorbar(1:n_par_it,corr_avg_sl,corr_CI_sl);
    grid on;
    set(gca,'XTick',1:n_par_it);
    set(gca,'XTickLabel',model.params.sigma_range);
    title('Variance explained per sigma ratio sl locked');
    xlabel('sigma ratio');
    ylabel('Variance explained');
    
    % Position ratio vs Variance explained for bb
    tmp_corr_bb = nan(n_par_it,size(chan_occ,1));
    for i=1:n_par_it
        tmp_corr_bb(i,:) = squeeze(all_corr(:,i,chan_occ,:,2));
    end
    corr_avg_bb = nanmean(tmp_corr_bb,2);
    corr_std_bb = nanstd(tmp_corr_bb,0,2);
    corr_ste_bb = corr_std_bb ./ sqrt(size(tmp_corr_bb,2));
    corr_CI_bb = 1.96 .* corr_ste_bb;
    fh_bb_VEparams = figure;
    plot(1:n_par_it,corr_avg_bb,'r','Linewidth',3);
    hold on; errorbar(1:n_par_it,corr_avg_bb,corr_CI_bb);
    grid on;
    set(gca,'XTick',1:n_par_it);
    set(gca,'XTickLabel',model.params.sigma_range);
    title('Variance explained per sigma ratio bb');
    xlabel('sigma ratio');
    ylabel('Variance explained');
    
    
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
    
    % Position ratio vs Variance explained for bb
    tmp_corr_bb = nan(n_par_it,size(chan_occ,1));
    for i=1:n_par_it
        tmp_corr_bb(i,:) = squeeze(all_corr(:,i,chan_occ,:,2));
    end
    corr_avg_bb = nanmean(tmp_corr_bb,2);
    corr_std_bb = nanstd(tmp_corr_bb,0,2);
    corr_ste_bb = corr_std_bb ./ sqrt(size(tmp_corr_bb,2));
    corr_CI_bb = 1.96 .* corr_ste_bb;
    fh_bb_VEparams = figure;
    plot(1:n_par_it,corr_avg_bb,'r','Linewidth',3);
    hold on; errorbar(1:n_par_it,corr_avg_bb,corr_CI_bb);
    grid on;
    set(gca,'XTick',1:n_par_it);
    set(gca,'XTickLabel',model.params.sigma_range);
    title('Variance explained per sigma ratio bb');
    xlabel('sigma ratio');
    ylabel('Variance explained');
    
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
    
    %     load(which('meg160_example_hdr.mat'))
    %     layout = ft_prepare_layout([],hdr);
    %     xpos = layout.pos(1:157,1);
    %     ypos = layout.pos(1:157,2);
    %     %chan_occ = find(ypos<1 & xpos<1);
    %     chan_occ = find(ypos<0 & xpos<1);
    %     % fh_allchan = mprfPlotHeadLayout(chan_occ);
    %     % saveas(fh_allchan,'Occchan_hp.jpg');
    
    tmp_corr = nan(n_par_it,size(chan_occ,1));
    for i=1:n_par_it
        tmp_corr(i,:) = squeeze(all_corr(1,i,chan_occ));
    end
    corr_avg = nanmean(tmp_corr,2);
    corr_std = nanstd(tmp_corr,0,2);
    corr_ste = corr_std ./ sqrt(size(tmp_corr,2));
    corr_CI = 1.96 .* corr_ste;
    corr_ci_occ = prctile(tmp_corr,[5 50 95],2);
    
    fh_sl_VEparams = figure;
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

save_dir = fullfile(main_dir, res_dir, rel_dir, ['Run_' type '_' opts.phs_metric '_' phfittype '_' cur_time]);
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

if ~isempty(fh_sl_VEparams)
    hgsave(fh_sl_VEparams,fullfile(save_dir,'Stimulus_locked_VEvsSig_occ'));
end

if ~isempty(fh_bb_VEparams)
    hgsave(fh_bb_VEparams,fullfile(save_dir,'Stimulus_locked_VEvsSig_occ'));
end

save(fullfile(save_dir, 'Results'),'results','model')

end
