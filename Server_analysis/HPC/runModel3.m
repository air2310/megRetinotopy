function runModel3()

    sub_sess_dir = pwd; 
    save_path = fullfile(sub_sess_dir, 'modeling', 'results');
    meg_data_file = fullfile(sub_sess_dir, 'data', 'meg', 'preproc', 'pp', 'epoched_data_hp_preproc_denoised.mat');    
    mprfSession_model_server(save_path,meg_data_file,sub_sess_dir,'3', 'n_cores', 4); 

end