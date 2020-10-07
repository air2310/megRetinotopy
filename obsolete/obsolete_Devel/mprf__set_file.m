function mprf__set_file(file_type, file_name)

global mprfSESSION
if isempty(mprfSESSION)
    error('mprfSESSION is empty')
end
   

switch lower(file_type)
    
    
    case {'vista_class', 'vista_class_file'}
        
        file_dir = mprf__get_directory(file_type);
        
        if strfind(file_type, '_file')
            
        else
            file_type = [file_type '_file'];
            
        end

        
        field_name = mprf__associate_field(file_type);
        eval(['mprfSESSION.' field_name '=''' fullfile(file_dir,file_name) ''';'])
        
        
        
    case {'vista_t1','vista_t1_file'}
        file_dir = mprf__get_directory(file_type);
        
        
        if strfind(file_type, '_file')
            
        else
            file_type = [file_type '_file'];
            
        end
        
        field_name = mprf__associate_field(file_type);
        eval(['mprfSESSION.' field_name '=''' fullfile(file_dir,file_name) ''';'])
        
    case {'rm_stim', 'rm_stim_file'}
        file_dir = mprf__get_directory('rm_stimulus');
        
        
        if strfind(file_type, '_file')
            
        else
            file_type = [file_type '_file'];
            
        end
        
        field_name = mprf__associate_field(file_type);
        eval(['mprfSESSION.' field_name '=''' fullfile(file_dir,file_name) ''';'])
        
        
    case 'prf_export_mat'
        file_dir = mprf__get_directory('prf_mat');
        
        
        field_name = mprf__associate_field(file_type);
        eval(['mprfSESSION.' field_name '=''' fullfile(file_dir,'exported_prf_params.mat') ''';'])
        
    case 'bs_anat'
        
        if strfind(file_type, '_file')
            
        else
            file_type = [file_type '_file'];
            
        end
        
        file_dir = mprf__get_directory(file_type);
        
        
        field_name = mprf__associate_field(file_type);
        eval(['mprfSESSION.' field_name '=''' fullfile(file_dir,file_name) ''';'])
        
        
    case 'fs_tag_to_idx'
        file_dir = mprf__get_directory('fs_surf_roi');
        field_name = mprf__associate_field(file_type);
        eval(['mprfSESSION.' field_name '=''' fullfile(file_dir,file_name) ''';'])
        
    case 'bs_tag_to_idx'
        
        file_dir = mprf__get_directory('bs_surf_roi');
        field_name =  mprf__associate_field(file_type);
        eval(['mprfSESSION.' field_name '=''' fullfile(file_dir,file_name) ''';'])
        
         
        
        
    otherwise
        fprint('Unknown file type: %s\n',file_type)
        
        
        
end





end