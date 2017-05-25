function mprfSession_prepare_session


load(fullfile(pwd,'mprfSESSION.mat'))
global mprfSESSION

main_dir = mprf__get_directory('main_dir');

if mprfSESSION.has.rm_model && ...
        mprfSESSION.has.vista_class && ...
        ~isempty(mprfSESSION.orig.paths.vista_path) && ...
        exist(mprfSESSION.orig.paths.vista_path,'dir');
    
    if isfield(mprfSESSION.has,'prf_exported') && exist(fullfile(mprf__get_directory('main_dir'),...
            mprf__get_file('prf_export_mat')),'file') && mprfSESSION.has.prf_exported
        
        fprintf('pRF parameters already imported and smoothed. Skipping.\n')

    else        
        
        if mprfSESSION.has.rm_stim_imported
            
            mprf__smooth_export_prf_parameters;
            
            
        else
            fprintf('We need the retinotopic model stimulus, but it has not been imported\n')
            fprintf('Importing it now\n')
            success = mprf__import_stimulus('rm_stim');
            
            if success
                mprf__smooth_export_prf_parameters;
                
            else
                fprintf('Error: could not import RM stimulus\n')
                return
                
                
            end
            
            
        end
    end
    
else
    fprintf(['Error: do not have the required data to export pRF parameters.\n'...
        'Please update session and import the required data'])
    return
    
    
end

success = mprf__export_roi_and_data_to_freesurfer_surface('both');

if success
    
    main_dir = mprf__get_directory('main_dir');
    
    bs_models = mprf__get_file('bs_headmodel');
    
    if isempty(bs_models.fname)
        fprint('Error: no Brainstorm head model found, import a head model and update this session\n')
        
    elseif length(bs_models.fname) == 1
        bs_surf_to_use = fullfile(bs_models.fpath,bs_models.fname(1).name);
        
    elseif length(bs_models.fname) > 1
        
        answer = listdlg('ListString',{bs_models.fname.name},...
            'PromptString','Multiple head models found, please select one',...
            'SelectionMode','Single');
        
        if isempty(answer)
            return
            
        end
        bs_surf_to_use = fullfile(bs_models.fpath, bs_models.fname(answer).name);
        
    end
    
    
    if mprfSESSION.has.bs_anat.subjectimage_T1 && exist(fullfile(main_dir, mprf__get_file('bs_anat')),'file')
        
        
    else
        fprintf('Error: no Brainstorm T1 found, please import Brainstorm anatomy folder and update this session\n')
        return
    end
    
    
    
    if success == 3 % Both ROI and data worked
        mprf__freesurfer_to_brainstorm('roi',bs_surf_to_use);
        mprf__freesurfer_to_brainstorm('prf',bs_surf_to_use);
        
    elseif success == 2 % Only data worked
        mprf__freesurfer_to_brainstorm('prf',bs_surf_to_use);
        
    elseif success == 1 % Only ROIs worked
        mprf__freesurfer_to_brainstorm('roi',bs_surf_to_use);
        
        
        
    else
        
    end
    
    
end


save(fullfile(main_dir, 'mprfSESSION.mat'),'mprfSESSION');
clear all














    
    
    






end