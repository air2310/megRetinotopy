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


%%


global mprfSESSION

if exist(fullfile(pwd,'mprfSESSION.mat'),'file')
    load('mprfSESSION.mat');
    
else
    error('Could not find mprfSESSION file. Please run from session folder')
end

main_dir = mprf__get_directory('main_dir');

[model, stimulus, syn_data] = mprfSession_model_gui;
[prf, bs, roi] = mprf__model_get_prf_params(model);

bs = mprf__get_lead_field(bs);
pred_resp = mprf__predicted_prf_response(model, stimulus, prf, roi);



cur_time = datestr(now);
cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';

save_dir = mprf__get_directory('model_predictions');
save(fullfile(main_dir, save_dir, ['model_predictions_' cur_time]));






