function mprf__load_mprfsession(action)

global mprfSESSION



switch lower(action)
    
    case 'init'
        
        cur_date = datestr(now);
        cur_date(cur_date == ' '| cur_date == ':' | cur_date == '-') = '_';
        
        mprfSESSION.init.date = cur_date;
        mprfSESSION.init.main_dir = pwd;
        mprfSESSION.root.root_path = mprfRootPath;


    case 'load'
        
        if ~exist(fullfile(pwd, 'mprfSESSION.mat'),'file')
            error('No mprfSESSION file found in current directory')
        end

        load(fullfile(pwd, 'mprfSESSION.mat'))
        
        if strcmpi(mprfSESSION.init.main_dir,pwd)
            
        else
            warning('Updating main directory to current directory')
            cur_date = datestr(now);
            cur_date(cur_date == ' '| cur_date == ':' | cur_date == '-') = '_';
            mprfSESSION.init.date = cur_date;
            mprfSESSION.init.main_dir = pwd;
            mprfSESSION.root.root_path = mprfRootPath;
        end


    otherwise
        error('Action not recognized')


end