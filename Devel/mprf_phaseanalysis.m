clear all;
% Averge phase from individual epochs of MEG data

%ref_ph_cri = 'mode';
%ref_ph_cri = 'amplitude';
ref_ph_cri = 'AmpPh';
refang_fit = 'h_varExp';
opts.metric = 'phaserefamplitude';
type = 'sl';
cur_time = mprf__get_cur_time;
main_dir = pwd;
save_dir = fullfile(main_dir, 'results', ['Run_' opts.metric '_' type '_'  ref_ph_cri '_' cur_time]);
%mkdir(save_dir);



%load('/home/akhi/Documents/Project/RetMEG/Subject_sessions/Akhil_tut/sub_004/MEG/data/preproc/epoched_data_hp_preproc_denoised.mat'); % Load MEG data (preprocessed)

load('data/meg/preproc/pp_run_R0774_RETMEG_Comb_5_12_17_on_26_Jul_2018_15_32_31/epoched_data_hp_preproc_denoised.mat');

% for sub 040
% data.data = epoched_filtered_preproc_denoised_data.data; % 1100 x 140 x 19 x 157
% (epoch x stim x repeat x sensor)
% for sub 004
 data.data = epoched_data.data;



%%
opts.stim_freq = 10;% Stimulus frequency
opts.samp_rate = 1000;% Sampling rate
Fs = opts.samp_rate;
% Get the size of the data:
sz = size(data.data);
opts.n_time = sz(1);
opts.n_bars = sz(2);
opts.n_reps = sz(3);
opts.n_chan = sz(4);

% Get the metric (defaults to amplitude)
%opts.metric = model.params.metric;
opts.idx = cell(1,1);

% Get the indices from the MEG epochs that correspond to the
% frequencies we need.
model.params.do_sl = 1;
model.params.do_bb = 0;
model.params.n_iterations = 1;
model.params.reliability_scans = 1;
model.params.reliability_split_half =1;
model.params.n_iterations_rel = 1;
model.params.n_cores = 4;
model.params.samp_rate = opts.samp_rate;
model.params.stim_freq = opts.stim_freq;
% 
% global mprfSESSION;
% mprfSESSION.init.main_dir = pwd;
% %model.type = 'reliability check';
 model.type = 'false';

% Define the periods (Eventually, this should come from somewhere else,
    % for example the stimulus files)
    periods.blank = [3:5 30:32 57:59 84:86 111:113 138:140]; % What frames where blanks?
    periods.blink = [1 2 28 29 55 56 82 83 109 110 136 137]; % During what frames could the subject blink (i.e. square in the center)
    periods.stim = setdiff(1:140,[periods.blink periods.blank]); % What where stimulus (i.e. bar position) frames?


%% (1)
% Check if the epochs are multiples of stimulus frequency (10 Hz)

%% (2)
% Select the most reliable sensor
if strcmpi(model.type,'reliability check')
    if model.params.reliability_scans || model.params.reliability_split_half
        
        if model.params.do_sl && model.params.do_bb
            fprintf('Running reliability for stimulus locked and broad band signal\n')
            data = [];
            data = mprf__do_reliability_analysis(model,'stimulus_locked',data);
            mprf__do_reliability_analysis(model,'broadband',data);
            
            
        elseif  model.params.do_sl && ~model.params.do_bb
            fprintf('Running reliability for stimiulus locked only\n')
            mprf__do_reliability_analysis(model,'stimulus_locked');
            
        elseif  ~model.params.do_sl && model.params.do_bb
            fprintf('Running reliability for broad band only\n')
            mprf__do_reliability_analysis(model,'broadband');
            I
        end
    end
end

% Results reliability check
% ********* Load from /home/akhi/Documents/Project/RetMEG/Subject_sessions/wl_subj_004_NEW/modeling/results/reliability_checks
[rel_fname, rel_fpath] = uigetfile('*.mat','Select reliability run results');
rel_results = load(fullfile(rel_fpath, rel_fname));

rel_mat = rel_results.results.split_half.corr_mat;
med_rel = nanmedian(rel_mat,2);
good_channels = med_rel > 0.1;

%% (3)
% Determine the phase and amplitude of the most reliable sensor

ft_data = mprf__fft_on_meg_data(data.data);% Fourier transform of the MEG data


if model.params.do_sl && model.params.do_bb
    [opts.idx{1}, opts.idx{2}] = mprf__get_freq_indices(true, true, opts);
    
elseif model.params.do_sl && ~model.params.do_bb
    opts.idx{1} = mprf__get_freq_indices(true, false, opts);
    
elseif model.params.do_bb && ~model.params.do_sl
    [~, opts.idx{1}] = mprf__get_freq_indices(false, true, opts);
    
else
    error('Unknown option')
    
end

%% Check if the epochs are multiples of stimulus frequency
% freq_del  = (1/1.1);
% freq_st = 0:freq_del:Fs/2;
% ft_amp = real(ft_data(:,58,5,103));
% 
% 
% tmp_orig_dt = data.data(1:100,58,5,103);
% figure, plot(tmp_orig_dt);
% 
% tmp_ft_data = ft_data(opts.idx{1},58,5,103);
% tmp_ft_amp = (2.*abs(tmp_ft_data))./opts.n_time;
% figure, plot(freq_st,tmp_ft_amp);

%%
% Determine the amplitude and VE of the fit to compute the most reliable
% sensor in terms of both high VE and high across scan reliability

opts.metric = 'amplitude';
[tseries_av_amp, tseries_std_amp, tseries_ste_amp, tseries_av_amp_rep] = mprf_computemetric(ft_data,opts,model);

[model_fname, model_fpath] = uigetfile('*.mat','Select model');
load(fullfile(model_fpath, model_fname));

% global mprfSESSION
% if exist(fullfile(pwd,'mprfSESSION.mat'),'file')
%     load('mprfSESSION.mat');
%     
% else
%     error('Could not find mprfSESSION file. Please run from session folder')
% end
% 
% if ~exist('read_curv','file')
%     tbUse({'vistasoft','retMEG'})
% end
% 
% main_dir = mprf__get_directory('main_dir');
% % Get the model options
% [model, stimulus, syn, channels] = mprfSession_model_gui;
% % Get the modeling parameters from the model variable:
% [prf, bs, roi] = mprf__model_get_prf_params(model);
% %load('model_predictions_02_Sep_2017_05_00_43.mat');
% iter = [];
% bs = mprf__get_lead_field(bs);
% pred_resp = mprf__predicted_prf_response(model, stimulus, prf, roi, iter);
% meg_resp = mprf__compute_predicted_meg_response(bs, pred_resp, channels);

% Fit the MEG predictions on the time series extracted above
n_it = model.params.n_iterations; % How many iterations do we need?
n_cores = model.params.n_cores; % How many cores are we going to use?
% Sizes (again...)
n_bars = size(meg_resp{1},1);
n_chan = size(meg_resp{1},2);
n_roi = size(meg_resp{1},3);
n_metric = size(tseries_av_amp,3);

preds_amp = nan(n_bars, n_chan, n_roi, n_metric);
mean_ve_amp = nan(n_chan, n_roi, n_metric);
fit_data_amp = nan(n_bars,n_chan,n_roi,n_metric);

for this_chan = 1:n_chan
    for this_roi = 1:n_roi
        for this_metric = 1:n_metric
            % Get MEG predictions:
            cur_pred = meg_resp{1}(:,this_chan);
            % Get current data:
            cur_data = tseries_av_amp(:,this_chan,this_metric);
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
            preds_amp(not_nan, this_chan, this_roi, this_metric) =  X * B;
            % Compute coefficient of determination (i.e. R square /
            % variance explained):
            cod_01 = 1- (var(cur_data(not_nan) - (X * B)) ./ var(cur_data(not_nan)));
            mean_ve_amp( this_chan, this_roi, this_metric) = cod_01;
            fit_data_amp(:, this_chan, this_roi, this_metric) = cur_data;
        end
    end
end

% Most reliable sensor based on the Amplitude and reliability
% good_channels = mean_ve_amp > 0.2;
% mst_rel_channel = find(med_rel == max(med_rel(good_channels)));
mst_rel_channel = find(mean_ve_amp == max(mean_ve_amp(good_channels)));

figPoint_2 = figure;
stem(mean_ve_amp,'k.','LineWidth',2,'MarkerSize',25);
hold on; stem(med_rel,'b.','LineWidth',2,'MarkerSize',25);
grid off;
title('VE and reliability');
xlabel('channels');
ylabel('VE');
legend('VE Amp','Reliability');
%saveas(figPoint_2,strcat(save_dir,'/Amp_vs_rel'),'fig');
%saveas(figPoint_2,strcat(save_dir,'/Amp_vs_rel'),'jpg');
%saveas(figPoint_2,strcat(save_dir,'/Amp_vs_rel'),'eps');


% sorting the reliability and Varaince explained based on the minimum value
tmp = [mean_ve_amp med_rel];
[B,sort_idx] = sort(min(tmp,[],2));
mean_ve_amp_sort = mean_ve_amp(sort_idx);
med_rel_sort = med_rel(sort_idx);
idx_nans = (isnan(mean_ve_amp_sort) | isnan(med_rel_sort));
sort_idx = sort_idx(~idx_nans);
rel_channels = sort_idx(end-10:end);
x_axes = 1:length(sort_idx(rel_channels));

figPoint_2_1 = figure;
stem(mean_ve_amp_sort,'k.','LineWidth',2,'MarkerSize',25);
hold on; stem(med_rel_sort,'b.','LineWidth',2,'MarkerSize',25);
set(gca,'XTick',[],'XTickLabel',[]);
grid off;
xlabel('channels');
ylabel('VE');
legend('VE Amp','Reliability','Location','NorthWest');
title('sorted VE and reliability');
xlabel('channels');
ylabel('VE');

%saveas(figPoint_2_1,strcat(save_dir,'/Amp_vs_rel_sorted'),'fig');
%saveas(figPoint_2_1,strcat(save_dir,'/Amp_vs_rel_sorted'),'jpg');
%saveas(figPoint_2_1,strcat(save_dir,'/Amp_vs_rel_sorted'),'eps');

figPoint_2_1_1 = figure(200);
stem(x_axes,mean_ve_amp(rel_channels),'filled','k.','LineWidth',2,'MarkerSize',25);
hold on; stem(x_axes,med_rel(rel_channels),'filled','b.','LineWidth',2,'MarkerSize',25);
set(gca,'XTick',x_axes,'XTickLabel',rel_channels);
grid on;
title('sorted VE and reliability (RELIABLE CHANNELS)');
xlabel('channels');
ylabel('VE');
legend('VE Amp','Reliability','Location','NorthWest');
%saveas(figPoint_2_1_1,strcat(save_dir,'/Amp_vs_rel_sorted_relchan'),'fig');
%saveas(figPoint_2_1_1,strcat(save_dir,'/Amp_vs_rel_sorted_relchan'),'jpg');
%saveas(figPoint_2_1_1,strcat(save_dir,'/Amp_vs_rel_sorted_relchan'),'eps');



%% Phase of the most reliable sensor
opts.metric = 'phase';
[tseries_av_ph, tseries_std_ph, tseries_ste_ph, tseries_av_ph_rep] = mprf_computemetric(ft_data,opts,model);

tseries_av_mst_rel_ph = tseries_av_ph(:,mst_rel_channel);
tseries_av_mst_rel_amp = tseries_av_amp(:,mst_rel_channel);


%% Reference phase computation per reliable channel


% Determine the reference phase
xbins = -3.14:0.314:3.14-0.314;
switch ref_ph_cri
    case {'mode'} % reference phase is the phase corresponding to the peak value of the distribution
        %mst_rel_ang = mode(tseries_av_mst_rel_ph(:));
        %[N,x] = hist(tseries_av_mst_rel_ph,20);
        xbins = -3.14:0.314:3.14-0.314;
        [N,x] = hist(tseries_av_mst_rel_ph,xbins);
        idx_max = find(N == max(N));
        mst_rel_ang = x(idx_max);
        PH_opt = mst_rel_ang;
        figPoint_3 = figure(100);
        rose(tseries_av_mst_rel_ph,xbins);
        hold on;
        x_refph = N(idx_max).*cos(xbins(idx_max));
        y_refph = N(idx_max).*sin(xbins(idx_max));
        plot([0 x_refph],[0 y_refph],'ro-');
        
        heading = strcat('polarplot ch',num2str(mst_rel_channel),{' '},num2str(rad2deg(PH_opt)));
        set(figPoint_3,'Name',heading{1});
        hold off;
        figPoint_3_1 = figure(101);
        hist(tseries_av_mst_rel_ph,xbins);
        % Polar histogram plot
        saveas(figPoint_3,strcat(save_dir,'/','PolarHist_phase_mstrelaiable'),'jpg');
        saveas(figPoint_3,strcat(save_dir,'/','PolarHist_phase_mstrelaiable'),'fig');
        saveas(figPoint_3,strcat(save_dir,'/','PolarHist_phase_mstrelaiable'),'eps');
        % Histogram plot
        saveas(figPoint_3_1,strcat(save_dir,'/','Hist_phase_mstrelaiable'),'jpg');
        saveas(figPoint_3_1,strcat(save_dir,'/','Hist_phase_mstrelaiable'),'fig');
        saveas(figPoint_3_1,strcat(save_dir,'/','Hist_phase_mstrelaiable'),'eps');
        
        
    case {'amplitude'} % ref phase is the phase corresponding to the maximum amplitude value
        max_amp_epoch = find(tseries_av_amp(:,mst_rel_channel) == max(tseries_av_amp(:,mst_rel_channel)));
        mst_rel_ang = tseries_av_ph(max_amp_epoch,mst_rel_channel);
        PH_opt = mst_rel_ang;
        figPoint_4 =figure;
        mprf_polarplot(tseries_av_mst_rel_amp, tseries_av_mst_rel_ph);
        hold on;
        x_refph = tseries_av_amp(max_amp_epoch,mst_rel_channel).*cos(mst_rel_ang);
        y_refph = tseries_av_amp(max_amp_epoch,mst_rel_channel).*sin(mst_rel_ang);
        plot([0 x_refph],[0 y_refph],'ro-');
        heading = strcat('polarplot ch',num2str(mst_rel_channel),{' '},num2str(rad2deg(PH_opt)));
        set(figPoint_4,'Name',heading{1});
        hold off;
        % Polar plot
        saveas(figPoint_4,strcat(save_dir,'/','Amplitude_vs_phase_mstrelaiable'),'jpg');
        saveas(figPoint_4,strcat(save_dir,'/','Amplitude_vs_phase_mstrelaiable'),'fig');
        saveas(figPoint_4,strcat(save_dir,'/','Amplitude_vs_phase_mstrelaiable'),'eps');
        
        
    case {'AmpPh'}
        
        tseries_av_rel_ph = tseries_av_ph(:,rel_channels);
        tseries_av_rel_amp = tseries_av_amp(:,rel_channels);
        ph_range = -3.14:0.0314:3.14;
        VE_all = nan(length(ph_range),size(tseries_av_rel_ph,2));
        ang_opt = nan(1,size(tseries_av_rel_ph,2));
        PH_opt = nan(1,size(tseries_av_rel_ph,2));
        VE_opt = nan(1,size(tseries_av_rel_ph,2));
        
        if strcmp(refang_fit,'h_var')
            for idx_ch = 1:size(tseries_av_rel_ph,2)
                tseries_av_rel_ph_ch =  tseries_av_rel_ph(:,idx_ch);
                tseries_av_rel_amp_ch =  tseries_av_rel_amp(:,idx_ch)./max(tseries_av_rel_amp(:,idx_ch));
                
                figPoint_2_2 = figure;
                %eval(strcat('figPoint_2_2_',num2str(idx_ch), '= figure'));
                mprf_polarplot(tseries_av_rel_amp_ch, tseries_av_rel_ph_ch);
                
                count_ang = 1;
                amp_diff_sq_sum = nan(size(ph_range));
                clear opt_diff_sq_sum VE_fit ref_ph_opt;
                for idx_ang = ph_range
                    ref_ph = idx_ang;
                    ph_diff = (tseries_av_rel_ph_ch - repmat(ref_ph,140,1));
                    amp_ref_ph = cos(ph_diff) .* tseries_av_rel_amp_ch;
                    
                    %amp_diff_sq_sum(:,count_ang) = nansum(amp_ref_ph.^2)./length(amp_ref_ph);
                    amp_diff_sq_sum(:,count_ang) = nanvar(amp_ref_ph);
                    
                    if notDefined('opt_diff_sq_sum')
                        opt_diff_sq_sum = amp_diff_sq_sum(:,count_ang);
                    end
                    if amp_diff_sq_sum(:,count_ang) >= opt_diff_sq_sum
                        ref_ph_opt = idx_ang;
                        VE_fit_opt = amp_diff_sq_sum(:,count_ang);
                        opt_diff_sq_sum = amp_diff_sq_sum(:,count_ang);
                        count_opt = count_ang;
                    end
                    count_ang = count_ang+1;
                end
                VE_all(:,idx_ch) = amp_diff_sq_sum';
                ang_opt(idx_ch) = count_opt;
                PH_opt(idx_ch) = ref_ph_opt;
                hold on;
                x = (amp_diff_sq_sum./max(amp_diff_sq_sum)).*cos(ph_range);
                y = (amp_diff_sq_sum./max(amp_diff_sq_sum)).*sin(ph_range);
                plot(x,y,'ko-','MarkerSize',1);
                hold on; plot([0 x(ang_opt(idx_ch))],[0 y(ang_opt(idx_ch))],'ro-');
                
                heading = strcat('polarplot ch',num2str(idx_ch),{' '},num2str(rad2deg(PH_opt(idx_ch))));
                set(figPoint_2_2,'Name',heading{1});
                hold off;
 %               saveas(figPoint_2_2,strcat(save_dir,'/',num2str(idx_ch),'polarplot_ch'),'fig');
 %               saveas(figPoint_2_2,strcat(save_dir,'/',num2str(idx_ch),'polarplot_ch',num2str(idx_ch)),'jpg');
 %               saveas(figPoint_2_2,strcat(save_dir,'/',num2str(idx_ch),'polarplot_ch',num2str(idx_ch)),'eps');
                
            end
            
        elseif strcmp(refang_fit,'h_varExp')
            VE_fit_alang = nan(size(ph_range,2),size(tseries_av_ph,2));
            num_chan_rel = size(tseries_av_ph,2);
            
            for idx_ch = 1:num_chan_rel
                %tseries_av_rel_ph_ch =  tseries_av_ph(:,idx_ch);
                %tseries_av_rel_amp_ch =  tseries_av_amp(:,idx_ch)./max(tseries_av_rel_amp(:,idx_ch));
                
                %figPoint_2_2 = figure;
                %eval(strcat('figPoint_2_2_',num2str(idx_ch), '= figure'));
                %mprf_polarplot(tseries_av_rel_amp_ch, tseries_av_rel_ph_ch);
                %pause();
                
                cur_ch = idx_ch;
                count_ang = 1;
                clear opt_ve VE_fit_opt;
                
                for idx_ang = ph_range
                    ref_ph = idx_ang;
                    opts.metric = 'phase ref amplitude';
                    %opts.metric = 'amplitude';
                    [tseries_av, tseries_std, tseries_ste] = mprf_computemetric(ft_data(:,:,:,cur_ch),opts,model,ref_ph);
                
                    % Compute the signed prediction
                    % check where the prediction is, remove the abs
                    
                    % Fit the MEG predictions on the time series extracted above
                    n_it = model.params.n_iterations; % How many iterations do we need?
                    n_cores = model.params.n_cores; % How many cores are we going to use?
                    % Sizes (again...)
                    n_bars = size(meg_resp{1},1);
                    n_chan = size(meg_resp{1},2);
                    n_roi = size(meg_resp{1},3);
                    n_metric = size(tseries_av,3);
                    
                    fit_data_opt = nan(n_bars,num_chan_rel);
                    preds_opt =  nan(n_bars,num_chan_rel);
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
                    
                    
                    VE_fit_alang(count_ang,idx_ch) = mean_ve(idx_ch);
                    count_ang = count_ang+1;
                end
                
                ang_opt(idx_ch) = count_opt;
                PH_opt(idx_ch) = ref_ph_opt;
                VE_opt(idx_ch) = VE_fit_opt;
                %hold on;
                
                x = (VE_fit_opt./max(VE_fit_opt)).*cos(ph_range);
                y = (VE_fit_opt./max(VE_fit_opt)).*sin(ph_range);
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
end

%% Plot of variance explained; MEG head plot - Amplitude and phase ref ampltude

% Variance explained
% For all unsorted channels
figPoint_ve_opt = figure;
stem(VE_opt,'filled','r.','LineWidth',2,'MarkerSize',25);
hold on;
stem(mean_ve_amp,'filled','k.','LineWidth',2,'MarkerSize',25);
xlabel('channels');
ylabel('VE');
legend('phase ref Amp','Amp','Location','NorthWest');
title('VE opt - Amp vs ph ref Amp (UNSORTED)');
saveas(figPoint_ve_opt,strcat(save_dir,'/','VE opt - Amp vs ph ref Amp'),'fig');
saveas(figPoint_ve_opt,strcat(save_dir,'/','VE opt - Amp vs ph ref Amp'),'jpg');
saveas(figPoint_ve_opt,strcat(save_dir,'/','VE opt - Amp vs ph ref Amp'),'eps');


% For all channels sorted
figPoint_ve_opt_sort = figure;
stem(VE_opt(sort_idx),'r.','LineWidth',2,'MarkerSize',25);
hold on;
stem(mean_ve_amp(sort_idx),'k.','LineWidth',2,'MarkerSize',25);
xlabel('channels');
ylabel('VE');
legend('phase ref Amp','Amp','Location','NorthWest');
title('VE opt - Amp vs ph ref Amp (SORTED)');
set(gca,'XTick',[],'XTickLabel',[]);
saveas(figPoint_ve_opt_sort,strcat(save_dir,'/','VE opt - Amp vs ph ref Amp sorted'),'fig');
saveas(figPoint_ve_opt_sort,strcat(save_dir,'/','VE opt - Amp vs ph ref Amp sorted'),'jpg');
saveas(figPoint_ve_opt_sort,strcat(save_dir,'/','VE opt - Amp vs ph ref Amp sorted'),'eps');

% For reliable channels
figPoint_ve_opt_sort_rel = figure;
stem(VE_opt(rel_channels),'filled','r.','LineWidth',2,'MarkerSize',25);
hold on;
stem(mean_ve_amp(rel_channels),'filled','k.','LineWidth',2,'MarkerSize',25);
xlabel('channels');
ylabel('VE');
legend('phase ref Amp','Amp','Location','NorthWest');
title('VE opt rel - Amp vs ph ref Amp');
x_axes = 1:length(rel_channels);
set(gca,'XTick',x_axes,'XTickLabel',rel_channels);
saveas(figPoint_ve_opt_sort_rel,strcat(save_dir,'/','VE_opt_Amp_vs_ph_ref_Amp_sorted_mstreliable'),'fig');
saveas(figPoint_ve_opt_sort_rel,strcat(save_dir,'/','VE_opt_Amp_vs_ph_ref_Amp_sorted_mstreliable'),'jpg');
saveas(figPoint_ve_opt_sort_rel,strcat(save_dir,'/','VE_opt_Amp_vs_ph_ref_Amp_sorted_mstreliable'),'eps');


% Head plot map
figPoint_rlch_map = mprfPlotHeadLayout(rel_channels);
saveas(figPoint_rlch_map,strcat(save_dir,'/','rel chan map'),'fig');
saveas(figPoint_rlch_map,strcat(save_dir,'/','rel chan map'),'jpg');
saveas(figPoint_rlch_map,strcat(save_dir,'/','rel chan map'),'eps');

HiVE_chan_amp = mean_ve_amp>0.2;
rel_channels_VE_amp = find(HiVE_chan_amp);
figPoint_rlchve_amp_map = mprfPlotHeadLayout(rel_channels_VE_amp);
saveas(figPoint_rlchve_amp_map,strcat(save_dir,'/','most VE amp map'),'fig');
saveas(figPoint_rlchve_amp_map,strcat(save_dir,'/','most VE amp map'),'jpg');
saveas(figPoint_rlchve_amp_map,strcat(save_dir,'/','most VE amp map'),'eps');

HiVE_chan = VE_opt>0.2;
rel_channels_VE = find(HiVE_chan);
figPoint_rlchve_map = mprfPlotHeadLayout(rel_channels_VE);
saveas(figPoint_rlchve_map,strcat(save_dir,'/','most VE map'),'fig');
saveas(figPoint_rlchve_map,strcat(save_dir,'/','most VE map'),'jpg');
saveas(figPoint_rlchve_map,strcat(save_dir,'/','most VE map'),'eps');

figPoint_phrf_map = figure;
megPlotMap(VE_opt,[0 0.6],figPoint_phrf_map,'jet','VE phase ref fit',[],[],'interpmethod','nearest');
saveas(figPoint_phrf_map,strcat(save_dir,'/','VE ph ref map'),'fig');
saveas(figPoint_phrf_map,strcat(save_dir,'/','VE ph ref map'),'jpg');
saveas(figPoint_phrf_map,strcat(save_dir,'/','VE ph ref map'),'eps');

figPoint_ampf_map = figure;
megPlotMap(mean_ve_amp,[0 0.6],figPoint_ampf_map,'jet','VE amp fit',[],[],'interpmethod','nearest');
saveas(figPoint_ampf_map,strcat(save_dir,'/','VE amp ref map'),'fig');
saveas(figPoint_ampf_map,strcat(save_dir,'/','VE amp ref map'),'jpg');
saveas(figPoint_ampf_map,strcat(save_dir,'/','VE amp ref map'),'eps');

figPoint_ampf_map = figure;
megPlotMap((VE_opt' - mean_ve_amp),[-0.5 0.6],figPoint_ampf_map,'jet','VE diff fit');
saveas(figPoint_ampf_map,strcat(save_dir,'/','VE diff ref map'),'fig');
saveas(figPoint_ampf_map,strcat(save_dir,'/','VE diff ref map'),'jpg');
saveas(figPoint_ampf_map,strcat(save_dir,'/','VE diff ref map'),'eps');

% Polar  histogram of phase values
figPoint_polhst = figure;
mprf_polarplot(ones(size(PH_opt)),PH_opt);
x_bl = ones(size(PH_opt)).* cos(PH_opt);
y_bl = ones(size(PH_opt)).* sin(PH_opt);
for ct = 1:size(PH_opt,2)
    hold on; plot([0 x_bl(ct)],[0 y_bl(ct)],'ro-');
end
saveas(figPoint_polhst,strcat(save_dir,'/','PolarHist_phase'),'jpg');
saveas(figPoint_polhst,strcat(save_dir,'/','PolarHist_phase'),'fig');
saveas(figPoint_polhst,strcat(save_dir,'/','PolarHist_phase'),'eps');


figPoint_polhst = figure;
mprf_polarplot(VE_opt,PH_opt);
saveas(figPoint_polhst,strcat(save_dir,'/','PolarHist_phase_VE'),'jpg');
saveas(figPoint_polhst,strcat(save_dir,'/','PolarHist_phase_VE'),'fig');
saveas(figPoint_polhst,strcat(save_dir,'/','PolarHist_phase_VE'),'eps');

% iio = (tseries_av_mst_rel_ph == 0);
% tseries_av_mst_rel_no0_ph = tseries_av_mst_rel_ph(~iio);
% mst_rel_ang = mode(tseries_av_mst_rel_no0_ph);
% hist(tseries_av_mst_rel);
% saveas(figPoint_3,'Hist_phase_mstrelaiable','jpg');

%% Phase referenced amplitude / signed predictions
for count_rel_ph = 1:length(PH_opt)
    
    mst_rel_ang = PH_opt(count_rel_ph);
    % Compute the phase referenced amplitude
    opts.metric = 'phaserefamplitude';
    %opts.metric = 'amplitude';
    [tseries_av, tseries_std, tseries_ste] = mprf_computemetric(ft_data,opts,model,mst_rel_ang);
    
    % Compute the signed prediction
    % check where the prediction is, remove the abs
    
    % Fit the MEG predictions on the time series extracted above
    n_it = model.params.n_iterations; % How many iterations do we need?
    n_cores = model.params.n_cores; % How many cores are we going to use?
    % Sizes (again...)
    n_bars = size(meg_resp{1},1);
    n_chan = size(meg_resp{1},2);
    n_roi = size(meg_resp{1},3);
    n_metric = size(tseries_av,3);
    
    preds = nan(n_bars, n_chan, n_roi, n_metric);
    mean_ve = nan(n_chan, n_roi, n_metric);
    fit_data = nan(n_bars,n_chan,n_roi,n_metric);
    
    for this_chan = 1:n_chan
        for this_roi = 1:n_roi
            for this_metric = 1:n_metric
                % Get MEG predictions:
                cur_pred = meg_resp{1}(:,this_chan);
                % Get current data:
                cur_data = tseries_av(:,this_chan,this_metric);
                % Exclude NANs
                not_nan = ~isnan(cur_pred(:)) & ~isnan(cur_data(:));
                
                % Make X matrix (predictor)
                if strcmpi(opts.metric, 'amplitude')
                    X = [ones(size(cur_pred(not_nan))) abs(cur_pred(not_nan))];
                elseif strcmpi(opts.metric,'phaserefamplitude')
                    X = [ones(size(cur_pred(not_nan))) (cur_pred(not_nan))];
                end
                % Compute Beta's:
                B = X \ cur_data(not_nan);
                % Store the predicted times series:
                preds(not_nan, this_chan, this_roi, this_metric) =  X * B;
                % Compute coefficient of determination (i.e. R square /
                % variance explained):
                cod_01 = 1- (var(cur_data(not_nan) - (X * B)) ./ var(cur_data(not_nan)));
                mean_ve( this_chan, this_roi, this_metric) = cod_01;
                fit_data(:, this_chan, this_roi, this_metric) = cur_data;
            end
        end
    end
    
    
    HiVE_chan = mean_ve > 0.2;
    sum(HiVE_chan);
    
    
    %% plots  
    %% Predictions and Original time series
    % Highest Variance explained channel for Phase ref amplitude
    ch_hve = find(mean_ve == max(mean_ve));
    figPoint_1 = figure;
    plot(tseries_av(:,ch_hve),'k','LineWidth',1);
    hold on; plot(preds(:,ch_hve),'r','LineWidth',2);
    grid on;
    xlabel('Time points')
    ylabel('A.U.')
    title(strcat('highest Variance explained CH ',':',num2str(ch_hve),{' '},'VE',':',num2str(max(mean_ve))));
    saveas(figPoint_1,strcat(save_dir,'/',num2str(count_rel_ph),'Highest VE'),'fig');
    saveas(figPoint_1,strcat(save_dir,'/',num2str(count_rel_ph),'Highest VE'),'jpg');
    saveas(figPoint_1,strcat(save_dir,'/',num2str(count_rel_ph),'Highest VE'),'eps');
    
    
    % Highest Variance explained channel for Amplitude
    figPoint_2 = figure;
    plot(tseries_av(:,mst_rel_channel),'k','LineWidth',1);
    hold on; plot(preds(:,mst_rel_channel),'r','LineWidth',2);
    grid on;
    xlabel('Time points')
    ylabel('A.U.')
    title(strcat('Most reliable CH',':',num2str(mst_rel_channel),{' '},'VE',':',num2str(mean_ve(mst_rel_channel))));
    saveas(figPoint_2,strcat(save_dir,'/',num2str(count_rel_ph),'Most reliable'),'fig');
    saveas(figPoint_2,strcat(save_dir,'/',num2str(count_rel_ph),'Most reliable'),'jpg');
    saveas(figPoint_2,strcat(save_dir,'/',num2str(count_rel_ph),'Most reliable'),'eps');
    
    
    % Most reliable channel for amplitude
    figPoint_5 = figure;
    plot(tseries_av_amp(:,ch_hve),'k','LineWidth',1);
    hold on; plot(preds_amp(:,ch_hve),'r','LineWidth',2);
    grid on;
    xlabel('Time points')
    ylabel('A.U.')
    title(strcat('highest Variance explained CH',':',num2str(ch_hve),{' '},'VE',':',num2str(mean_ve_amp(ch_hve))));
    saveas(figPoint_5,strcat(save_dir,'/',num2str(count_rel_ph),'Highest VE amp'),'fig');
    saveas(figPoint_5,strcat(save_dir,'/',num2str(count_rel_ph),'Highest VE amp'),'jpg');
    saveas(figPoint_5,strcat(save_dir,'/',num2str(count_rel_ph),'Highest VE amp'),'eps');
    
    % Most reliable channel for Phase ref amplitude
    figPoint_6 = figure;
    plot(tseries_av_amp(:,mst_rel_channel),'k','LineWidth',1);
    hold on; plot(preds_amp(:,mst_rel_channel),'r','LineWidth',2);
    grid on;
%     txt1 = strcat('VE =',num2str(mean_ve_amp(mst_rel_channel)));
%     txt2 = strcat('Ch = ',num2str(mst_rel_channel));
%     text(70,(2.7).*abs(tseries_av(70,ch_hve)),txt1);
%     text(70,2.5.*abs(tseries_av(70,ch_hve)),txt2);
    xlabel('Time points')
    ylabel('A.U.')
    title(strcat('Most reliable CH',':',num2str(mst_rel_channel),{' '},'VE',':',num2str(mean_ve_amp(mst_rel_channel))));
    saveas(figPoint_6,strcat(save_dir,'/',num2str(count_rel_ph),'Most reliable amp'),'fig');
    saveas(figPoint_6,strcat(save_dir,'/',num2str(count_rel_ph),'Most reliable amp'),'jpg');
    saveas(figPoint_6,strcat(save_dir,'/',num2str(count_rel_ph),'Most reliable amp'),'eps');
    
    
    % same channel fits
    chan_cur = rel_channels(count_rel_ph);
    figPoint_100 = figure;
    plot(tseries_av(:,chan_cur),'k','LineWidth',1);
    hold on; plot(preds(:,chan_cur),'r','LineWidth',2);
    grid on;
    xlabel('Time points')
    ylabel('A.U.')
    title(strcat('Current CH ',':',num2str(chan_cur),{' '},'VE',':',num2str(mean_ve(chan_cur))));
    saveas(figPoint_100,strcat(save_dir,'/',num2str(count_rel_ph),'current chan'),'fig');
    saveas(figPoint_100,strcat(save_dir,'/',num2str(count_rel_ph),'current chan'),'jpg');
    saveas(figPoint_100,strcat(save_dir,'/',num2str(count_rel_ph),'current chan'),'eps');
    
    % same channel for Amplitude
    figPoint_200 = figure;
    plot(tseries_av_amp(:,chan_cur),'k','LineWidth',1);
    hold on; plot(preds_amp(:,chan_cur),'r','LineWidth',2);
    grid on;
    xlabel('Time points')
    ylabel('A.U.')
    title(strcat('Current CH',':',num2str(chan_cur),{' '},'VE',':',num2str(mean_ve_amp(chan_cur))));
    saveas(figPoint_200,strcat(save_dir,'/',num2str(count_rel_ph),'current chan amp'),'fig');
    saveas(figPoint_200,strcat(save_dir,'/',num2str(count_rel_ph),'current chan amp'),'jpg');
    saveas(figPoint_200,strcat(save_dir,'/',num2str(count_rel_ph),'current chan amp'),'eps');
    
    %% Head plot - VE of individual channels on the head
    % Position of the sensors
    
    figPoint_rlch_map = mprfPlotHeadLayout(rel_channels);
    saveas(figPoint_rlch_map,strcat(save_dir,'/',num2str(count_rel_ph),'rel chan map'),'fig');
    saveas(figPoint_rlch_map,strcat(save_dir,'/',num2str(count_rel_ph),'rel chan map'),'jpg');
    saveas(figPoint_rlch_map,strcat(save_dir,'/',num2str(count_rel_ph),'rel chan map'),'eps');
    
    rel_channels_VE = find(HiVE_chan);
    figPoint_rlchve_map = mprfPlotHeadLayout(rel_channels_VE);
    saveas(figPoint_rlchve_map,strcat(save_dir,'/',num2str(count_rel_ph),'most VE map'),'fig');
    saveas(figPoint_rlchve_map,strcat(save_dir,'/',num2str(count_rel_ph),'most VE map'),'jpg');
    saveas(figPoint_rlchve_map,strcat(save_dir,'/',num2str(count_rel_ph),'most VE map'),'eps');
    
    figPoint_phrf_map = figure;
    megPlotMap(mean_ve,[0 0.6],figPoint_phrf_map,'jet','VE phase ref fit');
    saveas(figPoint_phrf_map,strcat(save_dir,'/',num2str(count_rel_ph),'VE ph ref map'),'fig');
    saveas(figPoint_phrf_map,strcat(save_dir,'/',num2str(count_rel_ph),'VE ph ref map'),'jpg');
    saveas(figPoint_phrf_map,strcat(save_dir,'/',num2str(count_rel_ph),'VE ph ref map'),'eps');    
    
    figPoint_ampf_map = figure;
    megPlotMap(mean_ve_amp,[0 0.6],figPoint_ampf_map,'jet','VE amp fit');
    saveas(figPoint_ampf_map,strcat(save_dir,'/',num2str(count_rel_ph),'VE amp map'),'fig');
    saveas(figPoint_ampf_map,strcat(save_dir,'/',num2str(count_rel_ph),'VE amp ref map'),'jpg');
    saveas(figPoint_ampf_map,strcat(save_dir,'/',num2str(count_rel_ph),'VE amp ref map'),'eps');
    
    %% Variance explained for the two fits
    % For all unsorted channels
    figPoint_8_1 = figure;
    stem(mean_ve,'filled','r.','LineWidth',2,'MarkerSize',25);
    hold on;
    stem(mean_ve_amp,'filled','k.','LineWidth',2,'MarkerSize',25);
    xlabel('channels');
    ylabel('VE');
    legend('phase ref Amp','Amp','Location','NorthWest');
    title('VE - Amp vs ph ref Amp (UNSORTED)');
    saveas(figPoint_8_1,strcat(save_dir,'/',num2str(count_rel_ph),'VE - Amp vs ph ref Amp'),'fig');
    saveas(figPoint_8_1,strcat(save_dir,'/',num2str(count_rel_ph),'VE - Amp vs ph ref Amp'),'jpg');
    saveas(figPoint_8_1,strcat(save_dir,'/',num2str(count_rel_ph),'VE - Amp vs ph ref Amp'),'eps');
    
    
    % For all channels sorted
    figPoint_8_2 = figure;
    stem(mean_ve(sort_idx),'r.','LineWidth',2,'MarkerSize',25);
    hold on;
    stem(mean_ve_amp(sort_idx),'k.','LineWidth',2,'MarkerSize',25);
    xlabel('channels');
    ylabel('VE');
    legend('phase ref Amp','Amp','Location','NorthWest');
    title('VE - Amp vs ph ref Amp (SORTED)');
    set(gca,'XTick',[],'XTickLabel',[]);    
    saveas(figPoint_8_2,strcat(save_dir,'/',num2str(count_rel_ph),'VE - Amp vs ph ref Amp sorted'),'fig');
    saveas(figPoint_8_2,strcat(save_dir,'/',num2str(count_rel_ph),'VE - Amp vs ph ref Amp sorted'),'jpg');
    saveas(figPoint_8_2,strcat(save_dir,'/',num2str(count_rel_ph),'VE - Amp vs ph ref Amp sorted'),'eps');
    
    % For reliable channels
    figPoint_8_3 = figure;
    stem(mean_ve(rel_channels),'filled','r.','LineWidth',2,'MarkerSize',25);
    hold on;
    stem(mean_ve_amp(rel_channels),'filled','k.','LineWidth',2,'MarkerSize',25);
    xlabel('channels');
    ylabel('VE');
    legend('phase ref Amp','Amp','Location','NorthWest');
    title('VE - Amp vs ph ref Amp');
    x_axes = 1:length(rel_channels);
    set(gca,'XTick',x_axes,'XTickLabel',rel_channels);
    saveas(figPoint_8_3,strcat(save_dir,'/',num2str(count_rel_ph),'VE_Amp_vs_ph_ref_Amp_sorted_mstreliable'),'fig');
    saveas(figPoint_8_3,strcat(save_dir,'/',num2str(count_rel_ph),'VE_Amp_vs_ph_ref_Amp_sorted_mstreliable'),'jpg');
    saveas(figPoint_8_3,strcat(save_dir,'/',num2str(count_rel_ph),'VE_Amp_vs_ph_ref_Amp_sorted_mstreliable'),'eps');
    
    
    results.mst_rel_chan = mst_rel_channel;
    results.mst_rel_ang =  mst_rel_ang;
    
    results.amp.mean_ve = mean_ve_amp;
    results.amp.fit_data = tseries_av_amp;
    results.amp.preds = preds_amp;
    
    results.ph.data = tseries_av_ph;
    
    results.ph_ref_amp.mean_ve = mean_ve;
    results.ph_ref_amp.fit_data = tseries_av;
    results.ph_ref_amp.preds = preds;
    
    save(strcat(save_dir,'/',num2str(count_rel_ph),'results'),'results');
    
end

%% Only for the figures

if strcmp(ref_ph_cri,'AmpPh')
    %ch = rel_channels(length(PH_opt));
    ch_all = [36,15,10,13,15,62,73,78];
    for count_chan = 1:length(ch_all)
        ch = ch_all(count_chan);
        fprintf('%d.',ch);
        fprintf('/n');
        for count_rel_ph = 1:length(PH_opt)
            
            mst_rel_ang = PH_opt(count_rel_ph);
            % Compute the phase referenced amplitude
            opts.metric = 'phaserefamplitude';
            %opts.metric = 'amplitude';
            [tseries_av, tseries_std, tseries_ste] = mprf_computemetric(ft_data,opts,model,mst_rel_ang);
            
            % Compute the signed prediction
            % check where the prediction is, remove the abs
            
            % Fit the MEG predictions on the time series extracted above
            n_it = model.params.n_iterations; % How many iterations do we need?
            n_cores = model.params.n_cores; % How many cores are we going to use?
            % Sizes (again...)
            n_bars = size(meg_resp{1},1);
            n_chan = size(meg_resp{1},2);
            n_roi = size(meg_resp{1},3);
            n_metric = size(tseries_av,3);
            
            preds = nan(n_bars, n_chan, n_roi, n_metric);
            mean_ve = nan(n_chan, n_roi, n_metric);
            fit_data = nan(n_bars,n_chan,n_roi,n_metric);
            
            for this_chan = 1:n_chan
                for this_roi = 1:n_roi
                    for this_metric = 1:n_metric
                        % Get MEG predictions:
                        cur_pred = meg_resp{1}(:,this_chan);
                        % Get current data:
                        cur_data = tseries_av(:,this_chan,this_metric);
                        % Exclude NANs
                        not_nan = ~isnan(cur_pred(:)) & ~isnan(cur_data(:));
                        
                        % Make X matrix (predictor)
                        if strcmpi(opts.metric, 'amplitude')
                            X = [ones(size(cur_pred(not_nan))) abs(cur_pred(not_nan))];
                        elseif strcmpi(opts.metric,'phaserefamplitude')
                            X = [ones(size(cur_pred(not_nan))) (cur_pred(not_nan))];
                        end
                        % Compute Beta's:
                        B = X \ cur_data(not_nan);
                        % Store the predicted times series:
                        preds(not_nan, this_chan, this_roi, this_metric) =  X * B;
                        % Compute coefficient of determination (i.e. R square /
                        % variance explained):
                        cod_01 = 1- (var(cur_data(not_nan) - (X * B)) ./ var(cur_data(not_nan)));
                        mean_ve( this_chan, this_roi, this_metric) = cod_01;
                        fit_data(:, this_chan, this_roi, this_metric) = cur_data;
                    end
                end
            end
            
            
            figPoint_1 = figure;
            plot(tseries_av(:,ch),'k','LineWidth',1);
            hold on; plot(preds(:,ch),'r','LineWidth',2);
            grid on;
            xlabel('Time points')
            ylabel('A.U.')
            title(strcat('CH ',':',num2str(ch),{' '},'VE',':',num2str(mean_ve(ch))));
            saveas(figPoint_1,strcat(save_dir,'/',num2str(count_rel_ph),'_',num2str(ch),'phrefAmp'),'fig');
            saveas(figPoint_1,strcat(save_dir,'/',num2str(count_rel_ph),'_',num2str(ch),'phrefAmp'),'jpg');
            saveas(figPoint_1,strcat(save_dir,'/',num2str(count_rel_ph),'_',num2str(ch),'phrefAmp'),'eps');
            
            % Most reliable channel for Phase ref amplitude
            figPoint_6 = figure;
            plot(tseries_av_amp(:,ch),'k','LineWidth',1);
            hold on; plot(preds_amp(:,ch),'r','LineWidth',2);
            grid on;
            xlabel('Time points')
            ylabel('A.U.')
            title(strcat('CH',':',num2str(ch),{' '},'VE',':',num2str(mean_ve_amp(ch))));
            saveas(figPoint_6,strcat(save_dir,'/',num2str(count_rel_ph),'_',num2str(ch),'Amp'),'fig');
            saveas(figPoint_6,strcat(save_dir,'/',num2str(count_rel_ph),'_',num2str(ch),'Amp'),'jpg');
            saveas(figPoint_6,strcat(save_dir,'/',num2str(count_rel_ph),'_',num2str(ch),'Amp'),'eps');
            
        end
    end
end