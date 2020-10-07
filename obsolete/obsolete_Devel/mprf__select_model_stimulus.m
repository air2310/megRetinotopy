function [stim_path, stim_file, stim_ext] = mprf__select_model_stimulus(stim_type)

stim_path = '';
stim_file = '';
stim_ext = '';

global FULL_STIM_PATH
main_dir = mprf__get_directory('main_dir');
make_stimulus = true;


stim_files = dir(fullfile(main_dir, mprf__get_directory(stim_type),'*.mat'));

stim_file_names = {stim_files.name};
stim_file_names = [stim_file_names 'Make new stimulus'];


answer = listdlg('ListString',stim_file_names,...
    'PromptString','Select existing stimulus, or make new one');

stim_opt = stim_file_names{answer};

if strcmpi(stim_opt,'make new stimulus')
    
else
    [stim_path, stim_file, stim_ext] = fileparts(fullfile(main_dir, mprf__get_directory(stim_type),stim_file_names{answer}));
    make_stimulus = false;
end


if make_stimulus && strcmpi(stim_type, 'rm_stimulus')
    
    errordlg('Option not implemented for retinotopic stimulus...')
    
elseif make_stimulus && strcmpi(stim_type, 'meg_imported_stim')
    
    answer = questdlg(sprintf('Do you want to use the retinotopic stimulus or an other stimulus?'),...
        'Select stimulus',...
        'Use retinotopic stimulus',...
        'Other stimulus',...
        'Other stimulus');
    
    if strcmpi(answer, 'use retinotopic stimulus')
        
        success = mprf__import_stimulus('rm_stim_as_pred');
        
        if success
            [stim_path, stim_file, stim_ext] = fileparts(FULL_STIM_PATH);
            clear FULL_STIM_PATH
            
        else
            
            
        end
        
    elseif strcmpi(answer,'other stimulus')
        
        
        success = mprf__import_stimulus('new_pred_stim');
        
        if success
            [stim_path, stim_file, stim_ext] = fileparts(FULL_STIM_PATH);
            clear FULL_STIM_PATH
            
        else
            
            
        end
        
    end
    
end


end













