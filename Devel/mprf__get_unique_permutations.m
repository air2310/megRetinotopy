function list = mprf__get_unique_permutations(samp_max, n_vals, max_perm)

init_vals = prod(samp_max - [0:n_vals-1]) ./ factorial(n_vals);
if init_vals > max_perm
    init_vals = max_perm;
end

list = nan(init_vals, n_vals);

n_unique = 0;
while n_unique ~= init_vals
    
    tmp = randperm(samp_max);
    cur_perm = sort(tmp(1:n_vals));
    is_on_list = false;
    if n_unique > 0
        for n = 1:n_unique
            
            if all(list(n,:) == cur_perm)
                is_on_list = true;
                break
                
            end
        end
        if ~is_on_list
            n_unique = n_unique +1;
            list(n_unique,:) = cur_perm;
            
        end
        
        
    else
        n_unique = n_unique +1;
        list(n_unique,:) = cur_perm;
        
    end
end

    
end







