function runModel1()

    addpath(genpath('/scratch/ek99/megRetinotopy'));
    sub_sess_dir = pwd; 
    save_path = fullfile(sub_sess_dir, 'modeling', 'results'); 
    meg_data_file = fullfile(sub_sess_dir, 'data', 'meg', 'raw'); 
    mprfSession_model_server(save_path,meg_data_file,sub_sess_dir,1); 
    exit()

end