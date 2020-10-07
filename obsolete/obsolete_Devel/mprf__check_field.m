function [result, no_fields] =  mprf__check_field(main_var, fields) %#ok<INUSL>
% Checks if a field exist and if it is empty...
% if field does not exist, return -1
% if field does exist, but is empty, return 0
% if field does exist, and is not empty, return 1


all_fields = strsplit(fields,'.');
no_fields = [];

for n = 1:length(all_fields)
    
    if n == 1;
        add_field = '';
        var_str = ['main_var' add_field];
        field_str = all_fields{n};
    else
        add_field = [add_field '.' all_fields{n-1}]; %#ok<AGROW>
        var_str = ['main_var' add_field];
        field_str = all_fields{n};
        
    end
    
    eval(['field_exist = isfield(' var_str ',''' field_str ''');'])
    if ~field_exist
        result = -1;
        no_fields = all_fields(n:end);
        return
    end
    
    if n == length(all_fields)
        add_field = [add_field '.' all_fields{n}]; %#ok<AGROW>
        var_str = ['main_var' add_field];
        
        eval(['field_empty = isempty(' var_str ');'])
        
        if field_empty
            result = 0;
        elseif ~field_empty
            result = 1;
        end
        
    end
    
end






end