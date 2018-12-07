function mprf_model = mprf_generate_model_params(model_type,varargin)

if nargin > 1
    addArg = varargin;
    if numel(addArg) == 1
        addArg=addArg{1};
    end
else
    addArg = [];
end

switch model_type
    case 'original (phase ref amplitude) (model fit)'
        mprf_model.type = 'Run Original model';
        mprf_model.cur_comment = 'Generate predictions based on the estimated pRF parameters, using the selected prediction stimulus';
        mprf_model.params = struct();
        mprf_model.params.sigma = 'sigma_smoothed'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
        mprf_model.params.x0 = 'x_smoothed';
        mprf_model.params.y0 = 'y_smoothed';
        mprf_model.params.beta = 'recomp_beta'; % Recomputed beta, max response, beta etc etc... fixed beta
        
        mprf_model.params.sigma_fix = NaN; % Absolute pRF size
        mprf_model.params.sigma_range = []; % Range of prf sizes, either absolute or proportional
        
        mprf_model.params.x0_fix = NaN;
        mprf_model.params.x0_range = [];
        
        mprf_model.params.y0_fix = NaN;
        mprf_model.params.y0_range = [];
        
        mprf_model.params.beta_fix = NaN;
        mprf_model.params.beta_range = [];
        
        mprf_model.params.roi_specific = 0; % 1 or 0
        mprf_model.params.fit_roi_specific = 0; % 1 or 0
        
        mprf_model.params.do_bb = 0;  % 1 or 0
        mprf_model.params.do_sl = 1;% 1 or 0
        mprf_model.params.stim_freq = 10; %Value
        
        mprf_model.params.reliability_split_half = 0; % 1 or 0; this is more like an analyses on the time series...
        mprf_model.params.reliability_scans = 0; % 1 or 0; this is also more like an analyses done on the time series, not really a modeling feature
        
        mprf_model.params.n_iterations_rel = NaN; % How many iterations for reliability check above?
        
        mprf_model.params.metric = 'phase ref amplitude'; % What metric to use? Amplitude, coherence, Stim/NBour freq etc ...
        mprf_model.params.metric_include_stim_freq = 0; % 1 of 0, include stimulus frequency when computing metric, i.e. Stim/Nbour freq
        mprf_model.params.metric_nbour_range = ''; % Value, i.e. the range of neighbouring frequencies to include
        
        mprf_model.params.fit_cross_validate = 0; % Do a cross validation when fitting
        mprf_model.params.fit_n_iterations_cv = ''; % how many iterations for this cross validation
        
        mprf_model.params.comb_lr_rois = 1;
        mprf_model.params.comb_dv_rois = 1;
        mprf_model.params.roi_mask = 1;
        
        mprf_model.params.ve_thr = 1;
        mprf_model.params.beta_thr = 1;
        mprf_model.params.ve_thr_vals = [0.1 inf];
        mprf_model.params.beta_thr_vals = [0 95];
        
        mprf_model.params.beta_equal_pred = 0;
        mprf_model.params.beta_equal_beta = 1;
        
        mprf_model.params.n_iterations = 1;
        mprf_model.params.n_cores = 1;
        mprf_model.params.samp_rate = 1000;
        
        mprf_model.params.phase_fit = 'model_fit';
        mprf_model.params.phase_fit_loo = 0;
        mprf_model.params.loo_type = 'Amp';
        
        
    case 'original (phase ref amplitude) (model fit) (leave one out)'
        mprf_model.type = 'Run Original model';
        mprf_model.cur_comment = 'Generate predictions based on the estimated pRF parameters, using the selected prediction stimulus';
        mprf_model.params = struct();
        mprf_model.params.sigma = 'sigma_smoothed'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
        mprf_model.params.x0 = 'x_smoothed';
        mprf_model.params.y0 = 'y_smoothed';
        mprf_model.params.beta = 'recomp_beta'; % Recomputed beta, max response, beta etc etc... fixed beta
        
        mprf_model.params.sigma_fix = NaN; % Absolute pRF size
        mprf_model.params.sigma_range = []; % Range of prf sizes, either absolute or proportional
        
        mprf_model.params.x0_fix = NaN;
        mprf_model.params.x0_range = [];
        
        mprf_model.params.y0_fix = NaN;
        mprf_model.params.y0_range = [];
        
        mprf_model.params.beta_fix = NaN;
        mprf_model.params.beta_range = [];
        
        mprf_model.params.roi_specific = 0; % 1 or 0
        mprf_model.params.fit_roi_specific = 0; % 1 or 0
        
        mprf_model.params.do_bb = 0;  % 1 or 0
        mprf_model.params.do_sl = 1;% 1 or 0
        mprf_model.params.stim_freq = 10; %Value
        
        mprf_model.params.reliability_split_half = 0; % 1 or 0; this is more like an analyses on the time series...
        mprf_model.params.reliability_scans = 0; % 1 or 0; this is also more like an analyses done on the time series, not really a modeling feature
        
        mprf_model.params.n_iterations_rel = NaN; % How many iterations for reliability check above?
        
        mprf_model.params.metric = 'phase ref amplitude'; % What metric to use? Amplitude, coherence, Stim/NBour freq etc ...
        mprf_model.params.metric_include_stim_freq = 0; % 1 of 0, include stimulus frequency when computing metric, i.e. Stim/Nbour freq
        mprf_model.params.metric_nbour_range = ''; % Value, i.e. the range of neighbouring frequencies to include
        
        mprf_model.params.fit_cross_validate = 0; % Do a cross validation when fitting
        mprf_model.params.fit_n_iterations_cv = ''; % how many iterations for this cross validation
        
        mprf_model.params.comb_lr_rois = 1;
        mprf_model.params.comb_dv_rois = 1;
        mprf_model.params.roi_mask = 1;
        
        mprf_model.params.ve_thr = 1;
        mprf_model.params.beta_thr = 1;
        mprf_model.params.ve_thr_vals = [0.1 inf];
        mprf_model.params.beta_thr_vals = [0 95];
        
        mprf_model.params.beta_equal_pred = 0;
        mprf_model.params.beta_equal_beta = 1;
        
        mprf_model.params.n_iterations = 1;
        mprf_model.params.n_cores = 1;
        mprf_model.params.samp_rate = 1000;
        
        mprf_model.params.phase_fit = 'model_fit';
        mprf_model.params.phase_fit_loo = 1;
        mprf_model.params.loo_type = 'Amp';
        
    case 'pRF size range (phase ref amplitude) (model fit) (leave one out)'
        mprf_model.type = 'pRF size range';
        mprf_model.cur_comment = 'Scaling the parameters gradually using a leave one out approach';
        mprf_model.params = struct();
        mprf_model.params.sigma = 'sigma_range (proportion, smoothed)'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
        mprf_model.params.x0 = 'x_smoothed';
        mprf_model.params.y0 = 'y_smoothed';
        mprf_model.params.beta = 'recomp_beta'; % Recomputed beta, max response, beta etc etc... fixed beta
        
        mprf_model.params.sigma_fix = NaN; % Absolute pRF size
        mprf_model.params.sigma_range = unique(round(logspace(log10(0.2),log10(10),20),1)); % Range of prf sizes, either absolute or proportional
        
        mprf_model.params.x0_fix = NaN;
        mprf_model.params.x0_range = [];
        
        mprf_model.params.y0_fix = NaN;
        mprf_model.params.y0_range = [];
        
        mprf_model.params.beta_fix = NaN;
        mprf_model.params.beta_range = [];
        
        mprf_model.params.roi_specific = 0; % 1 or 0
        mprf_model.params.fit_roi_specific = 0; % 1 or 0
        
        mprf_model.params.do_bb = 0;  % 1 or 0
        mprf_model.params.do_sl = 1;% 1 or 0
        mprf_model.params.stim_freq = 10; %Value
        
        mprf_model.params.reliability_split_half = 0; % 1 or 0; this is more like an analyses on the time series...
        mprf_model.params.reliability_scans = 0; % 1 or 0; this is also more like an analyses done on the time series, not really a modeling feature
        
        mprf_model.params.n_iterations_rel = NaN; % How many iterations for reliability check above?
        
        mprf_model.params.metric = 'phase ref amplitude'; % What metric to use? Amplitude, coherence, Stim/NBour freq etc ...
        mprf_model.params.metric_include_stim_freq = 0; % 1 of 0, include stimulus frequency when computing metric, i.e. Stim/Nbour freq
        mprf_model.params.metric_nbour_range = ''; % Value, i.e. the range of neighbouring frequencies to include
        
        mprf_model.params.fit_cross_validate = 0; % Do a cross validation when fitting
        mprf_model.params.fit_n_iterations_cv = ''; % how many iterations for this cross validation
        
        mprf_model.params.comb_lr_rois = 1;
        mprf_model.params.comb_dv_rois = 1;
        mprf_model.params.roi_mask = 1;
        
        mprf_model.params.ve_thr = 1;
        mprf_model.params.beta_thr = 1;
        mprf_model.params.ve_thr_vals = [0.1 inf];
        mprf_model.params.beta_thr_vals = [0 95];
        
        mprf_model.params.beta_equal_pred = 0;
        mprf_model.params.beta_equal_beta = 1;
        
        mprf_model.params.n_iterations = 1;
        mprf_model.params.n_cores = 1;
        mprf_model.params.samp_rate = 1000;
        
        mprf_model.params.phase_fit = 'model_fit';
        mprf_model.params.phase_fit_loo = 1;
        mprf_model.params.loo_type = 'Amp';
        
    case 'pRF position range (phase ref amplitude) (model fit) (leave one out)'
        mprf_model.type = 'position (x,y) range';
        mprf_model.cur_comment = 'Scaling the parameters gradually using a leave one out approach';
        mprf_model.params = struct();
        mprf_model.params.sigma = 'sigma_smoothed'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
        mprf_model.params.x0 = 'x range (proportion, smoothed)';
        mprf_model.params.y0 = 'y range (proportion, smoothed)';
        mprf_model.params.beta = 'recomp_beta'; % Recomputed beta, max response, beta etc etc... fixed beta
        
        mprf_model.params.sigma_fix = NaN; % Absolute pRF size
        mprf_model.params.sigma_range = []; % Range of prf sizes, either absolute or proportional
        
        mprf_model.params.x0_fix = NaN;
        mprf_model.params.x0_range = -pi:(pi/4):pi;
        
        mprf_model.params.y0_fix = NaN;
        mprf_model.params.y0_range = -pi:(pi/4):pi;
        
        mprf_model.params.beta_fix = NaN;
        mprf_model.params.beta_range = [];
        
        mprf_model.params.roi_specific = 0; % 1 or 0
        mprf_model.params.fit_roi_specific = 0; % 1 or 0
        
        mprf_model.params.do_bb = 0;  % 1 or 0
        mprf_model.params.do_sl = 1;% 1 or 0
        mprf_model.params.stim_freq = 10; %Value
        
        mprf_model.params.reliability_split_half = 0; % 1 or 0; this is more like an analyses on the time series...
        mprf_model.params.reliability_scans = 0; % 1 or 0; this is also more like an analyses done on the time series, not really a modeling feature
        
        mprf_model.params.n_iterations_rel = NaN; % How many iterations for reliability check above?
        
        mprf_model.params.metric = 'phase ref amplitude'; % What metric to use? Amplitude, coherence, Stim/NBour freq etc ...
        mprf_model.params.metric_include_stim_freq = 0; % 1 of 0, include stimulus frequency when computing metric, i.e. Stim/Nbour freq
        mprf_model.params.metric_nbour_range = ''; % Value, i.e. the range of neighbouring frequencies to include
        
        mprf_model.params.fit_cross_validate = 0; % Do a cross validation when fitting
        mprf_model.params.fit_n_iterations_cv = ''; % how many iterations for this cross validation
        
        mprf_model.params.comb_lr_rois = 1;
        mprf_model.params.comb_dv_rois = 1;
        mprf_model.params.roi_mask = 1;
        
        mprf_model.params.ve_thr = 1;
        mprf_model.params.beta_thr = 1;
        mprf_model.params.ve_thr_vals = [0.1 inf];
        mprf_model.params.beta_thr_vals = [0 95];
        
        mprf_model.params.beta_equal_pred = 0;
        mprf_model.params.beta_equal_beta = 1;
        
        mprf_model.params.n_iterations = 1;
        mprf_model.params.n_cores = 1;
        mprf_model.params.samp_rate = 1000;
        
        mprf_model.params.phase_fit = 'model_fit';
        mprf_model.params.phase_fit_loo = 1;
        mprf_model.params.loo_type = 'Amp';
        
    case 'scrambled (phase ref amplitude) (model fit) (leave one out)'
        mprf_model.type = 'scramble pRF parameters';
        mprf_model.cur_comment = 'Scrambling the parameters using a leave one out approach';
        mprf_model.params = struct();
        mprf_model.params.sigma = 'sigma_smoothed (scrambled,grouped)'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
        mprf_model.params.x0 = 'x_smoothed (scrambled,grouped)';
        mprf_model.params.y0 = 'y_smoothed (scrambled,grouped)';
        mprf_model.params.beta = 'recomp_beta (scrambled,grouped)'; % Recomputed beta, max response, beta etc etc... fixed beta
        
        mprf_model.params.sigma_fix = NaN; % Absolute pRF size
        mprf_model.params.sigma_range = []; % Range of prf sizes, either absolute or proportional
        
        mprf_model.params.x0_fix = NaN;
        mprf_model.params.x0_range = [];
        
        mprf_model.params.y0_fix = NaN;
        mprf_model.params.y0_range = [];
        
        mprf_model.params.beta_fix = NaN;
        mprf_model.params.beta_range = [];
        
        mprf_model.params.roi_specific = 0; % 1 or 0
        mprf_model.params.fit_roi_specific = 0; % 1 or 0
        
        mprf_model.params.do_bb = 0;  % 1 or 0
        mprf_model.params.do_sl = 1;% 1 or 0
        mprf_model.params.stim_freq = 10; %Value
        
        mprf_model.params.reliability_split_half = 0; % 1 or 0; this is more like an analyses on the time series...
        mprf_model.params.reliability_scans = 0; % 1 or 0; this is also more like an analyses done on the time series, not really a modeling feature
        
        mprf_model.params.n_iterations_rel = NaN; % How many iterations for reliability check above?
        
        mprf_model.params.metric = 'phase ref amplitude'; % What metric to use? Amplitude, coherence, Stim/NBour freq etc ...
        mprf_model.params.metric_include_stim_freq = 0; % 1 of 0, include stimulus frequency when computing metric, i.e. Stim/Nbour freq
        mprf_model.params.metric_nbour_range = ''; % Value, i.e. the range of neighbouring frequencies to include
        
        mprf_model.params.fit_cross_validate = 0; % Do a cross validation when fitting
        mprf_model.params.fit_n_iterations_cv = ''; % how many iterations for this cross validation
        
        mprf_model.params.comb_lr_rois = 1;
        mprf_model.params.comb_dv_rois = 1;
        mprf_model.params.roi_mask = 1;
        
        mprf_model.params.ve_thr = 1;
        mprf_model.params.beta_thr = 1;
        mprf_model.params.ve_thr_vals = [0.1 inf];
        mprf_model.params.beta_thr_vals = [0 95];
        
        mprf_model.params.beta_equal_pred = 0;
        mprf_model.params.beta_equal_beta = 1;
        
        mprf_model.params.n_iterations = 1000;
        mprf_model.params.n_cores = 1;
        mprf_model.params.samp_rate = 1000;
        
        mprf_model.params.phase_fit = 'model_fit';
        mprf_model.params.phase_fit_loo = 1;
        mprf_model.params.loo_type = 'Amp';
        
       
       case 'reliability'
        mprf_model.type = 'Reliability check';
        mprf_model.cur_comment = 'Do the reliability check';
        mprf_model.params = struct();
        mprf_model.params.sigma = 'sigma_smoothed'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
        mprf_model.params.x0 = 'x_smoothed';
        mprf_model.params.y0 = 'y_smoothed';
        mprf_model.params.beta = 'recomp_beta'; % Recomputed beta, max response, beta etc etc... fixed beta
        
        mprf_model.params.sigma_fix = NaN; % Absolute pRF size
        mprf_model.params.sigma_range = []; % Range of prf sizes, either absolute or proportional
        
        mprf_model.params.x0_fix = NaN;
        mprf_model.params.x0_range = [];
        
        mprf_model.params.y0_fix = NaN;
        mprf_model.params.y0_range = [];
        
        mprf_model.params.beta_fix = NaN;
        mprf_model.params.beta_range = [];
        
        mprf_model.params.roi_specific = 0; % 1 or 0
        mprf_model.params.fit_roi_specific = 0; % 1 or 0
            
        mprf_model.params.do_bb = 0;  % 1 or 0
        mprf_model.params.do_sl = 1;% 1 or 0
        mprf_model.params.stim_freq = 10; %Value
        
        mprf_model.params.reliability_split_half = 1; % 1 or 0; this is more like an analyses on the time series...
        mprf_model.params.reliability_scans = 1000; % 1 or 0; this is also more like an analyses done on the time series, not really a modeling feature
        
        mprf_model.params.n_iterations_rel = NaN; % How many iterations for reliability check above?
        
        mprf_model.params.metric = 'phase ref amplitude'; % What metric to use? Amplitude, coherence, Stim/NBour freq etc ...
        mprf_model.params.metric_include_stim_freq = 0; % 1 of 0, include stimulus frequency when computing metric, i.e. Stim/Nbour freq
        mprf_model.params.metric_nbour_range = ''; % Value, i.e. the range of neighbouring frequencies to include
        
        mprf_model.params.fit_cross_validate = 0; % Do a cross validation when fitting
        mprf_model.params.fit_n_iterations_cv = ''; % how many iterations for this cross validation
        
        mprf_model.params.comb_lr_rois = 1;
        mprf_model.params.comb_dv_rois = 1;
        mprf_model.params.roi_mask = 1;
        
        mprf_model.params.ve_thr = 1;
        mprf_model.params.beta_thr = 1;
        mprf_model.params.ve_thr_vals = [0.1 inf];
        mprf_model.params.beta_thr_vals = [0 95];
        
        mprf_model.params.beta_equal_pred = 0;
        mprf_model.params.beta_equal_beta = 1;
        
        mprf_model.params.n_iterations = 1;
        mprf_model.params.n_cores = 1;
        mprf_model.params.samp_rate = 1000;
        
        mprf_model.params.phase_fit = 'model_fit';
        mprf_model.params.phase_fit_loo = 0;
        mprf_model.params.loo_type = 'Amp'; 
end

mprf_model = rmProcessVarargin(mprf_model,addArg);

end


% To change the values of certain models
function mprf_model = rmProcessVarargin(mprf_model,vararg)
if ~exist('vararg','var') || isempty(vararg), return; end

fprintf(1,'[%s]:Resetting parameter:',mfilename);

for n=1:2:numel(vararg)
    data = vararg{n+1};
    fprintf(1,' %s,',vararg{n});
    
    switch lower(vararg{n})
        case {'n_iteration_scrambled'}
            if ischar(data)
                data = str2double(data);
            end
            mprf_model.params.n_iterations = data;
            
        case {'phase_fit_type'}
            mprf_model.params.phase_fit = data;
            
        case {'phase_fit_loo'}
            if ischar(data)
                data = str2double(data);
            end
            mprf_model.params.phase_fit_loo = data;
            
        case {'n_cores'}
            if ischar(data)
                data = str2double(data);
            end
            mprf_model.params.n_cores = data;
            
    end
    
end
fprintf(1,'.\n');
return
end