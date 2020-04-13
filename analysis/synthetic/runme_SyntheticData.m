%% runme_SyntheticData.m
% Input - Save path:  Path to save the mprf_params (model parameters)
%         meg_data_file_path: Folder where the preprocessed MEG data is saved  
%         sub_sess_dir: subject session folder
%         model_type: available model types - 
%                     'original (phase ref amplitude) (model fit)';
%                     'original (phase ref amplitude) (model fit) (leave one out)';
%                     'pRF size range (phase ref amplitude) (model fit) (leave one out)'
%                     'pRF position range (phase ref amplitude) (model fit) (leave one out)'
%                     'scrambled (phase ref amplitude) (model fit) (leave one out)'
%         varargin options available - 
%                     'n_iterations_scrambled'  
%                     'n_cores' - to run parallel for loops  
%                     'phase_fit_type' - method to determine reference phase - default: model_fit, can be changed to 'data_fit'
%                     'phase_fit_loo' - 1 indicates leave one out; 0 - not leave one out 


save_path     = fullfile(mprfRootPath,'Synthetic', 'results');
meg_data_file = fullfile(mprfRootPath,'Synthetic', 'meg_data.mat');
sub_sess_dir  = fullfile(mprfRootPath,'Synthetic');


%% Original fit - phase fit with model
model_type    = '1';
mprfSession_model_server(save_path,meg_data_file,sub_sess_dir,model_type)
