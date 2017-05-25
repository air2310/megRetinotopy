function [out] = getAllStructFields(field_in,field_name,out)
% field in: the current field (just the entire struct)
% field name: name of the current field (name of the entire struct)
% out: current output variable (empty cell)



fn = fieldnames(field_in);
for n = 1:length(fn)
    cur_name = [field_name '.' fn{n}];
    if isstruct(field_in.(fn{n}))
        cur_field = field_in.(fn{n});
        out = getAllStructFields(cur_field,cur_name,out);
    else
        out = [out cur_name];        
    end
end