function out_file = mprf__get_file(file_type)

global mprfSESSION
if isempty(mprfSESSION)
    error('mprfSESSION is empty')
end

out_file = [];
switch lower(file_type)
    
    
    case {'vista_class', 'vista_class_file'}
        if strfind(file_type, '_file')
            
        else
            file_type = [file_type '_file'];
            
        end
        
        
        field_name = mprf__associate_field(file_type);
        eval(['out_file = mprfSESSION.' field_name ';'])
        
        
        
    case {'vista_t1','vista_t1_file'}
        if strfind(file_type, '_file')
            
        else
            file_type = [file_type '_file'];
            
        end
        
        field_name = mprf__associate_field(file_type);
        eval(['out_file = mprfSESSION.' field_name ';'])
        
        
    case {'rm_stim','rm_stim_file'}
        if strfind(file_type, '_file')
            
        else
            file_type = [file_type '_file'];
            
        end
        
        field_name = mprf__associate_field(file_type);
        eval(['out_file = mprfSESSION.' field_name ';'])
    
    case 'prf_export_mat'
        
        field_name = mprf__associate_field(file_type);
        eval(['out_file = mprfSESSION.' field_name ';'])
        
        
    case 'bs_headmodel'
        % Allow for multiple possible head models
        
        main_dir = mprf__get_directory('main_dir');
        file_dir = mprf__get_directory('bs_headmodel');
        
        fpath = fullfile(main_dir, file_dir);
        tmp = dir(fullfile(fpath,'*.mat'));
        
        out_file.fpath = fpath;
        out_file.fname = tmp;
        
        
    case 'bs_anat'
        if strfind(file_type, '_file')
            
        else
            file_type = [file_type '_file'];
            
        end
        
        
        field_name = mprf__associate_field(file_type);
        eval(['out_file = mprfSESSION.' field_name ';'])
        
        
    case 'fs_tag_to_idx'
        
        field_name = mprf__associate_field(file_type);
        eval(['out_file = mprfSESSION.' field_name ';'])
        
        
    case 'bs_tag_to_idx'
        
        field_name = mprf__associate_field(file_type);
        eval(['out_file = mprfSESSION.' field_name ';'])
        
        
        
    otherwise
        fprintf('Unknown file type: %s\n',file_type)
        
        
end









end