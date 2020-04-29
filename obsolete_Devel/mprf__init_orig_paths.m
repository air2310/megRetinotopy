function mprf__init_orig_paths

global mprfSESSION

if isempty(mprfSESSION)
    error('mprfSESSION is empty')
end

[root_field, fields_to_add] = mprf__get_orig_path_fields;

default_value = '';


for n = 1:length(fields_to_add)
   
   cur_field = [root_field '.' fields_to_add{n}];
   result = mprf__check_field(mprfSESSION, cur_field);
   
   if result == -1
       eval(['mprfSESSION.' root_field '.' fields_to_add{n} ' = ''' default_value ''';'])
       fprintf('Field created: %s \n', cur_field)
   elseif result == 0
       fprintf('Field already exists, but is empty: %s \n', cur_field)
       
   elseif result == 1
       fprintf('Field already exists and is not empty: %s\n', cur_field)
       
   end
    
    
end

end





















