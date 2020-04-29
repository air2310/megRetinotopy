


stim_type = 'rm_stimulus';
make_stimulus = true;
stim_files = dir(fullfile(handles.main_dir, mprf__get_directory(stim_type),'*.mat'));


stim_file_names = {stim_files.name};
stim_file_names = [stim_file_names 'Make new stimulus'];


answer = listdlg('ListString',stim_file_names,...
    'PromptString','Select existing stimulus, or make new one');

stim_opt = stim_file_names{answer};

if strcmpi(stim_opt,'make new stimulus')
    
else
    stim_file_to_use = fullfile(fdir,stim_file_names{answer});
    make_stimulus = false;
end


if make_stimulus && strcmpi(stim_type, 'rm_stimulus')
    
    errordlg('Option not implemented for retinotopic stimulus...')
    
elseif make_stimulus && strcmpi(stim_type, 'prd_stimulus')
    
    answer = questdlg(sprintf('Do you want to use the retinotopic stimulus or an other stimulus?'),...
        'Select stimulus',...
        'Use retinotopic stimulus',...
        'Other stimulus',...
        'Other stimulus');
    
    if strcmpi(answer, 'use retinotopic stimulus')
        
        success = mprf__import_stimulus('rm_stim_as_pred');
        
    elseif strcmpi(answer,'other stimulus')
        
        
        success = mprf__import_stimulus('new_pred_stim');
        
        
        [s_name, s_path] = uigetfile('*.mat','Please select stimulus file to use');
        cur_dir = pwd;
        cd(s_path);
        [g_name, g_path] = uigetfile('*.mat','Please select corresponding grid file');
        cd(cur_dir);
        
        stim_file = fullfile(s_path, s_name);
        grid_file = fullfile(g_path, g_name);
        
        stim_file_to_use = mprfSessionPrepareStimulus({stim_file, grid_file});
        
    end
    
end
















