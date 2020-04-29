function mprfDeleteExcept(varargin)


if all(cellfun(@ischar, varargin))
    
    w = whos;
    for ii = 1:length(w)
        
        switch w(ii).name
            
            case varargin
                
                % do nothing
                
            otherwise
                
                clear(w(ii).name)
                
        end
        
    end
    
    
else
    error('All inputs must be char');
    
    
end
end


