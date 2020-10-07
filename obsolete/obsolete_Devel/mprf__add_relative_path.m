function mprf__add_relative_path(data_type, rel_path, overwrite) %#ok<INUSD>

global mprfSESSION
if isempty(mprfSESSION) % We do not expect mprfSESSION to be empty here, right
   error('mprfSESSION is empty')
    
end

if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = false;
end

    field_name = mprf__associate_field(data_type);
    
    if overwrite
        eval(['mprfSESSION.' field_name ' = rel_path;'])
    else
        result = mprf__check_field(mprfSESSION,field_name);
        
        if result > 0
            warning('Fields %s already exist and are not empty', rel_path)
        end
        
        answer = questdlg(sprintf('Field %s is not empty. Overwrite?',rel_path),...
            'Field is not empty','Keep','Overwrite','Keep');
        
        if strcmpi(answer,'keep')
            return
        elseif strcmi(answer,'overwrite')
            eval(['mprfSESSION.' field_name ' = rel_path;'])
        end
        
        
    end
    
    



end