function all_fields = getAllFieldsFromStruct(main_var)
% Recurrent function that returns all the field names of a struct,
% regardless of their 'depth'.
all_fields = getAllStructFields(main_var,'',{});

end