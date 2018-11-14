function runModel4()

    sub_sess_dir = pwd; 
    save_path = fullfile(sub_sess_dir, 'modeling', 'results');
    d = dir(fullfile(sub_sess_dir, 'data', 'meg', 'preproc', '*', 'epoched_data_hp_preproc_denoised.mat'));    
    meg_data_file = fullfile(d.folder, d.name); 
    mprfSession_model_server(save_path,meg_data_file,sub_sess_dir,'4'); 

end