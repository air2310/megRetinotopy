function out = mprfUpdateSession(nargs_ppair, varargin)

global mprfSESSION

out = 0;

if ~isnumeric(nargs_ppair)
    error('First argument must be an integer (N arguments per pair')
end

if mod(length(varargin),nargs_ppair)
    error('Varargin must have %d arguments for each input',nargs_ppair)
end
    
cur_date = datestr(now);
cur_date(cur_date == ' '| cur_date == ':' | cur_date == '-') = '_';



for n = 1:nargs_ppair:length(varargin)
    
    
    
    
    switch lower(varargin{n})
        
        
        
        case 'init'
            
            mprfSESSION.init.date = cur_date;
            mprfSESSION.init.main_dir = pwd;

            
        case 'add_orig_paths'
            out = max([exist(varargin{n+2},'dir') exist(varargin{n+2},'file')]); 
            if out
                mprfSESSION.orig.(varargin{n+1}) = varargin{n+2};
            end

            
            
            
        case 'remove_orig_paths'
            
            
        case 'add_source_path'
            mprfSESSION.source.(varargin{n+1}) = varargin{n+2};
            
        case 'remove_source_path'
            
            
            
        case 'add_source_data'
            mprfSESSION.source.(varargin{n+1}) = varargin{n+2};
        
        case 'add_root_path'
            out = max([exist(varargin{n+2},'dir') exist(varargin{n+2},'file')]);
            if out
                mprfSESSION.root.(varargin{n+1}) = varargin{n+2};
                
            end
            
        case 'remove_source_data'
            
            
            
            
            
        otherwise
            
            warning('Unknown action %s. Skipping',varargin{n})
            
            
            
            
    end
end

end











