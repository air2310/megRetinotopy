%function mprfSession_run_model

% % %% TO DO:
% % - Add options to the GUI that allow to combine ROIs
% dorsal and ventral, left and right
% - Set the broad band range in an
% intelligent way -> look at Eline's project for defaults
%
%   Allow run model to set some additional thresholds, like VE and beta
%   threshold
%
% Make gui to load model files that shows the most important options
% describing the file
% Make gui to visualize predicted responses on mesh, also adjust
% mprfSessionRenderDataOnBrainstormSurface and other functions to be
% compatible with new code set up
% Adjust MEG data preprocessing script, make gui as well with some
% options??
%
% Create code that replicates the first attempt and expand on this when
% necessary


%%
global mprfSESSION
if exist(fullfile(pwd,'mprfSESSION.mat'),'file')
    load('mprfSESSION.mat');
    
else
    error('Could not find mprfSESSION file. Please run from session folder')
end

if ~exist('notDefined','file')
    tbUse({'vistasoft','retMEG'})
end

main_dir = mprf__get_directory('main_dir');

[model, stimulus, syn, channels] = mprfSession_model_gui;

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
    
elseif strcmpi(model.type,'scramble pRF parameters') || ...
        strcmpi(model.type,'run original model') || ...
        strcmpi(model.type,'fix prf size') || ...
        strcmpi(model.type, 'prf size range')
    
    [prf, bs, roi] = mprf__model_get_prf_params(model);

    
    n_range = sum([prf.sigma.type.range
        prf.x0.type.range
        prf.y0.type.range
        prf.beta.type.range]);
    
    n_scramble = sum([prf.sigma.type.scrambled
        prf.x0.type.scrambled
        prf.y0.type.scrambled
        prf.beta.type.scrambled]);
    
    iter = [];
    
    if n_range >= 1 && n_scramble >= 1
       error('Both a range manipulation and a scrambling is not implemented') 
       
    elseif n_range == 0 && n_scramble == 0
            
    elseif n_range == 1 && n_scramble == 0
        tmp_fields = fieldnames(prf);
        iter.method = 'range';
        for n = 1:length(tmp_fields)
            if prf.(tmp_fields{n}).type.range
                
                iter.var = tmp_fields{n};
                iter.n = size(prf.(tmp_fields{n}).val,2);
                
            end
        end
        
    elseif n_scramble >= 1 && n_range == 0;
       iter.method = 'scramble';        
    elseif n_range > 1
        error('Range not implemented for multiple variables')
        
    else
        error('Unknown option')
        
        
    end
    bs = mprf__get_lead_field(bs);
    
    % BASICALLY WE SHOULD DO SOMETHING DIFFERENT IF WE'RE SCRAMBLING HERE.
    % THESE FUNCTIONS ACTUALLY STORE THE PREDICTIONS, WHICH IS FINE IF
    % THERE ARE ONLY 15 OF THEM, NOT WHEN THERE ARE A THOUSAND...
    
    pred_resp = mprf__predicted_prf_response(model, stimulus, prf, roi, iter);
    meg_resp = mprf__compute_predicted_meg_response(bs, pred_resp, channels);
    
    cur_time = datestr(now);
    cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
    
    save_dir = mprf__get_directory('model_predictions');
    save(fullfile(main_dir, save_dir, ['model_predictions_' cur_time]),...
        'prf','bs','roi','model','stimulus','pred_resp','syn','meg_resp','channels');
    
    if syn.do
        out_name = mprf__make_synthetic_data_set(syn, meg_resp, cur_time,channels);
        
    end
    
    if syn.pp
        mprfSession_preprocess_meg_data(true,out_name);
        
    end
end