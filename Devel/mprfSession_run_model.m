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

main_dir = mprf__get_directory('main_dir');

[model, stimulus, syn, channels] = mprfSession_model_gui;

if strcmpi(model.type,'reliability check')
    if model.params.reliability_scans || model.params.reliability_split_half
        
        if model.params.do_sl
            mprf__do_reliability_analysis(model,'stimulus_locked');
        end
        
        if model.params.do_bb
            mprf__do_reliability_analysis(model,'broadband');
            
        end
    end
    
else
    
    [prf, bs, roi] = mprf__model_get_prf_params(model);
    
    bs = mprf__get_lead_field(bs);
    pred_resp = mprf__predicted_prf_response(model, stimulus, prf, roi);
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