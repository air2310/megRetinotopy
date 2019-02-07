function mprfSession_model_server(save_path,meg_data_file,sub_sess_dir,model_type,varargin)
% mprfSession_model_server -  Script to run the meg retinotopy model on the
% server
% Input - Save path:  Path to save the mprf_params (model parameters)
%         meg_data_file_path: Folder where the preprocessed MEG data is saved  
%         sub_sess_dir: subject session folder
%         model_type: available model types - 
%                     'original (phase ref amplitude) (model fit)';
%                     'original (phase ref amplitude) (model fit) (leave one out)';
%                     'pRF size range (phase ref amplitude) (model fit) (leave one out)'
%                     'pRF position range (phase ref amplitude) (model fit) (leave one out)'
%                     'scrambled (phase ref amplitude) (model fit) (leave one out)'
%				      'reliability'		
%         varargin options available - 
%                     'n_iteration_scrambled'  
%                     'n_cores' - to run parallel for loops  
%                     'parpool_sys' - to create parallel pool for specific cluster profile - add switch case to the open_parallel_pool.m file                    
%                     'phase_fit_type' - method to determine reference phase - default: model_fit, can be changed to 'data_fit'
%                     'phase_fit_loo' - 1 indicates leave one out; 0 - not leave one out 
%                     'roi_mask' - 1 - include ROI mask, 0 - no ROI mask 
%                     

allModels = {'original (phase ref amplitude) (model fit)', ...
                 'original (phase ref amplitude) (model fit) (leave one out)', ...
                 'pRF size range (phase ref amplitude) (model fit) (leave one out)', ...
                 'pRF position range (phase ref amplitude) (model fit) (leave one out)', ...
                 'scrambled (phase ref amplitude) (model fit) (leave one out)'...
                 'reliability'};

if ~any(strcmp(allModels, model_type))
    model_type = allModels{str2double(model_type)};
end

global mprfSESSION
mprfSESSION = make_mprfSession(sub_sess_dir);

% Model to run
if ~exist('model_type','var') || isempty(model_type)
    model_type = 'original (phase ref amplitude) (model fit) (leave one out)';
end


% Adding the variable arguments to model
if nargin > 3
    addArg = varargin;
else
    addArg = [];
end



% Creates a file containing all the information about the model to run
[params,~] = mprf_make_params(save_path,model_type,addArg); 


% Script to run the model based on the parameters from the params_file
mprfSession_run_model_server(params,meg_data_file,0);

end