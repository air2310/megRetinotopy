function def_params = mprf__get_default_model_parameters(model_type)


def_params = init__default_params_struct;

%%% ADD SCRAMBLE STIMULUS OPTION AND N SCRAMBLE ITERATIONS
switch lower(model_type)
    
    
    case 'run original model'
        
        
    case 'equal weight'
        
        def_params.beta = 'fixed recomp_beta (absolute)'; % Recomputed beta, max response, beta etc etc... fixed beta
        def_params.beta_fix = '1';
        
        
    case 'scramble prf parameters'
        
        %%% ADD SCRAMBLE ITERATIONS
        
        def_params.sigma = 'sigma_smoothed (scrambled)'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
        def_params.x0 = 'x_smoothed (scrambled)';
        def_params.y0 = 'y_smoothed (scrambled)';
        def_params.beta = 'recomp_beta (scrambled)'; % Recomputed beta, max response, beta etc etc... fixed beta
        def_params.n_iterations_scramble = '1000';
        
    case 'fix prf size'
        
        def_params.sigma = 'fixed sigma (absolute)'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
        def_params.sigma_fix = '0.5'; % Absolute pRF size
        
    case 'prf size range'
        
        def_params.sigma = 'sigma range (proportion)'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
        def_params.sigma_range = '0.2:0.2:3'; % Range of prf sizes, either absolute or proportional
        
        def_params.x0 = 'x_smoothed';
        def_params.y0 = 'y_smoothed';
        def_params.beta = 'recomp_beta';
        
    case 'reliability check'
        def_params.reliability_split_half = 1; % 1 or 0; this is more like an analyses on the time series...
        def_params.reliability_scans = 1; % 1 or 0; this is also more like an analyses done on the time series, not really a modeling feature
        def_params.n_iterations_rel = '1000'; % How many iterations for reliability check above?
        
    case 'fit separate roi predictions'
        
        def_params.beta = 'fixed recomp_beta (absolute)'; % Recomputed beta, max response, beta etc etc... fixed beta
        def_params.beta_fix = 1;
        
        def_params.roi_specific = 1; % 1 or 0
        def_params.fit_roi_specific = 1; % 1 or 0

        def_params.fit_cross_validate = 1;
        def_params.fit_n_iterations_cv = '1000';

    otherwise
        error('Unknown model type %s.\n',model_type)
        
        
end

end

function def_params = init__default_params_struct

def_params.sigma = 'sigma_smoothed'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
def_params.x0 = 'x_smoothed';
def_params.y0 = 'y_smoothed';
def_params.beta = 'recomp_beta'; % Recomputed beta, max response, beta etc etc... fixed beta

def_params.sigma_fix = ''; % Absolute pRF size
def_params.sigma_range = ''; % Range of prf sizes, either absolute or proportional

def_params.x0_fix = '';
def_params.x0_range = '';

def_params.y0_fix = '';
def_params.y0_range = '';

def_params.beta_fix = '';
def_params.beta_range = '';

def_params.roi_specific = 0; % 1 or 0
def_params.fit_roi_specific = 0; % 1 or 0

def_params.do_bb = 1;  % 1 or 0
def_params.do_sl = 1;% 1 or 0
def_params.stim_freq = '10'; %Value

def_params.reliability_split_half = 0; % 1 or 0; this is more like an analyses on the time series...
def_params.reliability_scans = 0; % 1 or 0; this is also more like an analyses done on the time series, not really a modeling feature

def_params.n_iterations_rel = ''; % How many iterations for reliability check above?

def_params.metric = 'amplitude'; % What metric to use? Amplitude, coherence, Stim/NBour freq etc ...
def_params.metric_include_stim_freq = 0; % 1 of 0, include stimulus frequency when computing metric, i.e. Stim/Nbour freq
def_params.metric_nbour_range = ''; % Value, i.e. the range of neighbouring frequencies to include

def_params.fit_cross_validate = 0; % Do a cross validation when fitting
def_params.fit_n_iterations_cv = ''; % how many iterations for this cross validation

def_params.comb_lr_rois = 1;
def_params.comb_dv_rois = 1;
def_params.roi_mask = 1;

def_params.ve_thr = 1;
def_params.beta_thr = 1;
def_params.ve_thr_vals = {'0.1' 'inf'};
def_params.beta_thr_vals = {'0' '95'};

def_params.beta_equal_pred = 0;
def_params.beta_equal_beta = 1;

def_params.n_iterations_scramble = '';
def_params.n_cores = '1';
def_params.samp_rate = '1000';

end
