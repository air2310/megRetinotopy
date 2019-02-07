function mprf_model = getModelParams(modelType,varargin)
% Function that loads requested model parameters

% NOTE THIS FUNCTION IS STILL UNDER CONSTRUCTIOIN -- DO NOT USE

% mprf_model = getModelParams(modelType,varargin)

% INPUTS
% modelType    : string with an integer (default = '1');
%                    '1' - standard / original model (comparing data with prediction)
%                    '2' - prf size range model (comparing data with model prediction for different prf sizes)
%                    '3' - prf position range model (comparing data with model prediction for different prf positions)
%                    '4' - splif half reliability check using standard / original model
%                    '5' - scramble prf parameters across vertices (keeping x,y,sigma together)
% varargin option :  to be added..

if nargin > 1
    addArg = varargin;
    if numel(addArg) == 1
        addArg=addArg{1};
    end
else
    addArg = [];
end

mprf_model = struct();

%% Default parameters

% pRF params
mprf_model.params.sigma         = 'sigma_smoothed'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
mprf_model.params.x0            = 'x_smoothed';
mprf_model.params.y0            = 'y_smoothed';
mprf_model.params.beta          = 'recomp_beta'; % Recomputed beta, max response, beta etc etc... fixed beta

% Fix or change pRF parameters
mprf_model.params.sigma_fix     = NaN; % Absolute pRF size
mprf_model.params.sigma_range   = [];  % Range of prf sizes, either absolute or proportional
mprf_model.params.x0_fix        = NaN; % use fmri defined x0
mprf_model.params.x0_range      = [];  % range x0 parameter by subtracting/adding what value?
mprf_model.params.y0_fix        = NaN; % use fmri defined y0
mprf_model.params.y0_range      = [];  % range y0 parameter by subtracting/adding what value?
mprf_model.params.beta_fix      = NaN; % use fmri defined beta
mprf_model.params.beta_range    = [];  % range beta parameter by subtracting/adding what value?

% ROIs
mprf_model.params.roi_specific     = 0; % Use specific roi for ?? either 1 or 0
mprf_model.params.fit_roi_specific = 0; % Use specific roi for model fit? either 1 or 0
mprf_model.params.comb_lr_rois     = 1; % Combine left and right ROIs? either 1 or 0
mprf_model.params.comb_dv_rois     = 1; % Combine dorsal and ventral ROIs? either 1 or 0
mprf_model.params.roi_mask         = 1; % Use roi mask instead of individual rois? either 1 or 0

% MEG Data options
mprf_model.params.do_bb                     = 0;     % Use broadband power as signal of interest? either 1 or 0
mprf_model.params.do_sl                     = 1;     % Use stimulus-locked signal as signal of interest? either 1 or 0
mprf_model.params.stim_freq                 = 10;    % [Hz] Integer to define frequency of stimulus locked peak
mprf_model.params.samp_rate                 = 1000;  % [Hz] Sample rate
mprf_model.params.phase_fit                 = 'model_fit'; % Use what data for computing the reference phase? can be 'model_fit' or 'data_fit'
mprf_model.params.phase_fit_loo             = 1;     % Use leave-one-out crossvalidaton method to establish reference phase? either 1 or 0
mprf_model.params.loo_type                  = 'Amp'; % What MEG signal to use to calculate reference phase with leave one out crossvalidation? 'Amp' or ???
mprf_model.params.metric                    = 'phase ref amplitude'; % What metric to use? Can be: Amplitude, coherence, Stim/NBour freq etc ...
mprf_model.params.metric_include_stim_freq  = 0;      % either 1 of 0, include stimulus frequency when computing metric, i.e. Stim/Nbour freq
mprf_model.params.metric_nbour_range        = '';     % Value, i.e. the range of neighbouring frequencies to include

% Compute split half reliability?
mprf_model.params.reliability_split_half    = 0;      % either 1 or 0; this is more like an analyses on the time series...
mprf_model.params.reliability_scans         = 0;      % either 1 or 0; this is also more like an analyses done on the time series, not really a modeling feature
mprf_model.params.n_iterations_rel          = NaN;    % How many iterations for reliability check above?

mprf_model.params.fit_cross_validate        = 0;      % Do a cross validation when fitting? either 1 or 0
mprf_model.params.fit_n_iterations_cv       = '';     % How many iterations for this cross validation?

% Thresholding options
mprf_model.params.ve_thr                    = 1;
mprf_model.params.beta_thr                  = 1;
mprf_model.params.ve_thr_vals               = [0.1 inf];
mprf_model.params.beta_thr_vals             = [0 95];

mprf_model.params.beta_equal_pred           = 0;
mprf_model.params.beta_equal_beta           = 1;

% Computation options
mprf_model.params.n_iterations              = 1;
mprf_model.params.n_cores                   = 1;


switch modelType
    case '1' % Original / default
        mprf_model.type = 'Original model';
        mprf_model.cur_comment = 'Generate predictions based on the estimated pRF parameters, using the selected prediction stimulus';
         
    case '2' % Change pRF size in prediction
        mprf_model.type = 'pRF size range';
        mprf_model.cur_comment = 'Scaling the parameters gradually using a leave one out approach';
        
        mprf_model.params.sigma = 'sigma_range (proportion, smoothed)'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
        mprf_model.params.sigma_range = unique(round(logspace(log10(0.2),log10(10),20),1)); % Range of prf sizes, either absolute or proportional

    case '3' % Change pRF position in prediction
        mprf_model.type = 'position (x,y) range';
        mprf_model.cur_comment = 'Scaling the parameters gradually using a leave one out approach';
      
        mprf_model.params.x0 = 'x range (proportion, smoothed)';
        mprf_model.params.y0 = 'y range (proportion, smoothed)';
        
        mprf_model.params.x0_range = -pi:(pi/4):pi;  
        mprf_model.params.y0_range = -pi:(pi/4):pi;
        
    case '4' % Scramble pRF parameters across vertices (keeping x,y, sigma together/grouped)
        mprf_model.type = 'scramble pRF parameters';
        mprf_model.cur_comment = 'Scrambling the parameters using a leave one out approach';
       
        mprf_model.params.sigma = 'sigma_smoothed (scrambled,grouped)'; % sigma smoothed, sigma, fixed sigma (absolute), fixed sigma (proportion), sigma range (absolute), sigma range (proportion)
        mprf_model.params.x0 = 'x_smoothed (scrambled,grouped)';
        mprf_model.params.y0 = 'y_smoothed (scrambled,grouped)';
        mprf_model.params.beta = 'recomp_beta (scrambled,grouped)'; % Recomputed beta, max response, beta etc etc... fixed beta 
        mprf_model.params.n_iterations = 1000;
       
       case '5' % Do a split half reliability check using the original/default model
        mprf_model.type = 'Reliability check';
        mprf_model.cur_comment = 'Do the reliability check';
     
        mprf_model.params.reliability_split_half = 1; % 1 or 0; this is more like an analyses on the time series...
        mprf_model.params.reliability_scans = 1000; % 1 or 0; this is also more like an analyses done on the time series, not really a modeling feature
        mprf_model.params.phase_fit_loo = 0;
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