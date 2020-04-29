function success = mprf__import_stimulus(w_stim)

% We need the mprfSESSION file:
global mprfSESSION
if isempty(mprfSESSION)
    error('mprfSESSION is empty')
end

if strcmpi(w_stim,'rm_stim_as_pred') || strcmpi(w_stim, 'new_pred_stim')
   global FULL_STIM_PATH 
    
end


success = false;

switch lower(w_stim)
    
    case 'rm_stim'
        
        if ~mprfSESSION.has.vista_path || ~mprfSESSION.has.rm_model
            fprintf('It appears that the necessary files to import the retinotopic stimulus have not been imported. Skipping\n')
            
            
        else
            
            
            [~, out_file] = fileparts(mprfSESSION.orig.paths.rm_model);
            out_name = ['rm_stim_' out_file '.mat'];
            
            dest_path = fullfile(mprf__get_directory('main_dir'), mprf__get_directory('rm_stimulus'), out_name);
            
            if exist(dest_path,'file')
                answer = questdlg(sprintf('%s already exists. Overwrite?',out_name));
                
                if strcmpi(answer,'no') || strcmpi(answer,'cancel')
                    return
                end
                
                
                
            end
            
            % This function relies on several mrVista functions, so navigate to the
            % subject's mrVista directory
            cd(mprfSESSION.orig.paths.vista_path);
            
            
            % Check if we still can get to the mprfSESSION code directory:
            if exist('mprfLoadRMStimulus_02','file')
                
            else
                addpath(mprfSESSION.root.root_path);
                if exist('mprfLoadRMStimulus_02','file')
                    
                else
                    error('Could not add root directory to path')
                    
                end
            end
            
            
            % We need the dataTYPES variable
            mrGlobals;
            
            % We need a volume view:
            hvol = initHiddenGray;
            
            % Try to work out with data type the retinotopic model was based on by
            % matching the directories of the pRF model path to the dataTYPES:
            dt = {dataTYPES.name};
            dt_match = false(size(dt));
            
            rm_parse = strsplit(mprfSESSION.orig.paths.rm_model,'/');
            
            
            for n = 1:length(dt);
                cur_dt = dt{n};
                
                if any(strcmpi(cur_dt,rm_parse));
                    dt_match(n) = true;
                end
            end
            
            % If only one match has been found, we can be confident that we found the
            % correct dataTYPE, otherwise, bring up a dialog for the user to select it.
            if sum(dt_match) == 1
                use_dt = dt{dt_match};
            else
                use_dt = listdlg('listdlg',dt,...
                    'SelectionMode','single',...
                    'PrompString','Please select the data type to which the retinotopic model belongs');
            end
            
            % Close gracefully when no valid data type is selected
            if isempty(use_dt)
                return
                
            end
            
            % Set the volume view to the current data type and add the RM model
            hvol = viewSet(hvol,'curdt',lower(use_dt));
            hvol = rmSelect(hvol,1,mprfSESSION.orig.paths.rm_model);
            
            rmparams = viewGet(hvol,'rmparams');
            
            rm_stim.im = rmparams.stim.images_unconvolved;
            rm_stim.im_conv = rmparams.analysis.allstimimages';
            
            rm_stim.window = rmparams.stim.stimwindow;
            rm_stim.X = rmparams.analysis.X;
            rm_stim.Y = rmparams.analysis.Y;
            
            [rm_stim.full_x, rm_stim.full_y] = meshgrid(unique(rm_stim.X), unique(rm_stim.Y));
            [~, rm_stim.full_im] = rmStimulusMatrix(rmparams,[],[],false,false);
            
            
            save(dest_path,'rm_stim')
            mrvCleanWorkspace;
            
            if exist(dest_path,'file')
                mprfSESSION.has.rm_stim_imported = true;
                fprintf('Imported %s\n',dest_path);
                
                
            else
                mprfSESSION.has.rm_stim_imported = false;
                fprintf('Coule not import RM stimulus\n');
                
            end
            
            mprf__set_file('rm_stim',out_name)
            
            cd(mprf__get_directory('main_dir'))
            
            success = true;
            
        end
    case 'meg_stim'
        if ~mprfSESSION.has.meg_stim
            
            fprintf('It appears that the necessary files to import the MEG stimulus have not been imported. Skipping\n')
            
            
        else
            out_name = 'meg_stimulus.mat';
            dest_path = fullfile(mprf__get_directory('main_dir'), ...
                mprf__get_directory('meg_imported_stim'),out_name);
            
            
            if exist(dest_path,'file')
                answer = questdlg(sprintf('%s already exists. Overwrite?',out_name));
                
                if strcmpi(answer,'no') || strcmpi(answer,'cancel')
                    return
                end
                
                
                
            end
            
            
            indir = dir(fullfile(mprf__get_directory('main_dir'), ...
                mprf__get_directory('meg_stim'),...
                '*.mat'));
            
            fnames = {indir.name};
            tmp  = cellfun(@(x) strfind(x,'retinotopy_stimulus'), ...
                lower(fnames),'UniformOutput',false);
            fidx = find(~cellfun(@isempty,tmp), 1, 'first');
            
            if isempty(fidx)
                [fname, fpath] = uigetfile('*.mat','Please select stimulus file to use as a template');
                
            else
                fname = fnames{fidx};
                fpath = fullfile(mprf__get_directory('main_dir'), ...
                    mprf__get_directory('meg_stim'));
                
            end
            
            templ_file = fullfile(fpath, fname);
            
            fprintf('Using %s as MEG stimulus file\n',fname)
            
            tmp  = cellfun(@(x) strfind(x,'grid'), ...
                lower(fnames),'UniformOutput',false);
            grid_idx = find(~cellfun(@isempty,tmp), 1, 'first');
            
            if isempty(grid_idx)
                [gname, gpath] = uigetfile('*.mat','Please select Grid file to use');
                
                answer = questdlg('Apparently, a grid file has not been imported (properly). Import file now?');
                
                if strcmpi(answer,'yes')
                    yes_to_all = false;
                    [yes_to_all, success] = mprf__import_data(fullfile(gpath, gname),'meg_grid','MEG_grid.mat',yes_to_all);
                    
                elseif strcmpi(answer,'no')
                    
                elseif strcmpi(answer,'cancel')
                    
                    
                end
                
                
            else
                gname = fnames{grid_idx};
                gpath = fullfile(mprf__get_directory('main_dir'), ...
                    mprf__get_directory('meg_stim'));
                
                
                
            end
            
            grid_file = fullfile(gpath, gname);
            fprintf('Using %s as Grid file\n',gname)
            
            load(templ_file);
            load(grid_file);
            
            x_range = [min(grid.Xd(:)) max(grid.Xd(:))];
            y_range = [min(grid.Yd(:)) max(grid.Yd(:))];
            
            stim_range = double([min(stimulus.images(:)) max(stimulus.images(:))]);
            tmp = stimulus.images(:,:,1);
            bk = double(mode(tmp(:)));
            
            new_period = stimulus.seq(stimulus.trigSeq > 0);
            
            im_out_size = [101 101];
            meg_stim.full_im = zeros([im_out_size, length(new_period)]);
            
            srate = max(grid.Xd(:)) ./ floor(max(im_out_size) / 2);
            
            for n = 1:length(new_period)
                tmp_im = double(stimulus.images(:,:,new_period(n)));
                tmp_im = imresize(tmp_im,im_out_size,'nearest');
                meg_stim.full_im(:,:,n) = double(ceil(abs(tmp_im - bk) ./ max(abs(stim_range - bk)))).*(srate.^2);
                
            end
            
            [meg_stim.full_x, meg_stim.full_y] = meshgrid(linspace(x_range(1), x_range(2), im_out_size(1)),...
                linspace(y_range(1), y_range(2), im_out_size(2)));
            
            keep_window = sum(meg_stim.full_im,3) > 0;
            
            meg_stim.im = reshape(meg_stim.full_im,[],size(meg_stim.full_im,3));
            meg_stim.im =  meg_stim.im(keep_window(:),:);
            meg_stim.X = meg_stim.full_x(:);
            meg_stim.Y = meg_stim.full_y(:);
            meg_stim.X = meg_stim.X(keep_window);
            meg_stim.Y = meg_stim.Y(keep_window);
            
            
            
            save(dest_path,'meg_stim');
            
            if exist(dest_path,'file')
                mprfSESSION.has.meg_stim_imported = true;
                fprintf('Imported %s\n',dest_path);
                success = true;
                
            else
                mprfSESSION.has.meg_stim_imported = false;
                fprintf('Could not import MEG stimulus\n');
                
            end
          
            
        end
        
         
    case 'rm_stim_as_pred'
        
        main_dir = mprf__get_directory('main_dir');
        
        if mprfSESSION.has.rm_stim_imported
            
            rm_files = dir(fullfile(main_dir, mprf__get_directory('rm_stimulus'),'*.mat'));
            
            cur_time = datestr(now);
            cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
            
            dest_name = ['rm_stim_as_pred_' cur_time '.mat'];
            
            
            if length(rm_files) >= 1
                
                if length(rm_files) > 1

                
                end
                
                source_file = fullfile(main_dir, mprf__get_directory('rm_stimulus'), rm_files(1).name);
                dest_file = fullfile(main_dir, mprf__get_directory('meg_imported_stim'), dest_name);
                
                str_cmnd = sprintf('cp %s %s',source_file, dest_file);
                system(str_cmnd)
                
                if exist(dest_file,'file')
                    success = true;
                    FULL_STIM_PATH = dest_file;
                else
                    
                    
                end
                
                
                
            elseif isempty(rm_files)
                success_import = mprf__import_stimulus('rm_stim');
                success_set_pred = mprf__import_stimulus('rm_stim_as_pred');
                
                if success_import && success_set_pred
                    success = true;
                    
                else
                    
                    
                end
                
                
            end
            
        else
            
            
            
        end
        
        
    case 'new_pred_stim'
        
        main_dir = mprf__get_directory('main_dir');
        
        
        [s_name, s_path] = uigetfile('*.mat','Please select stimulus file to use');
        cur_dir = pwd;
        cd(s_path);
        [g_name, g_path] = uigetfile('*.mat','Please select corresponding grid file');
        cd(cur_dir);
        
        stim_file = fullfile(s_path, s_name);
        grid_file = fullfile(g_path, g_name);
        
        cur_time = datestr(now);
        cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
        
        dest_name = ['pred_stim' cur_time '.mat'];
        dest_path = fullfile(main_dir, mprf__get_directory('meg_imported_stim'));
        dest_file = fullfile(main_dir, mprf__get_directory('meg_imported_stim'), dest_name);
        
        
        fprintf('Using %s as MEG stimulus file\n',s_name)
        
        fprintf('Using %s as Grid file\n',g_name)
        
        load(stim_file);
        load(grid_file);
        
        x_range = [min(grid.Xd(:)) max(grid.Xd(:))];
        y_range = [min(grid.Yd(:)) max(grid.Yd(:))];
        
        stim_range = double([min(stimulus.images(:)) max(stimulus.images(:))]);
        tmp = stimulus.images(:,:,1);
        bk = double(mode(tmp(:)));
        
        new_period = stimulus.seq(stimulus.trigSeq > 0);
        
        im_out_size = [101 101];
        meg_stim.full_im = zeros([im_out_size, length(new_period)]);
        
        srate = max(grid.Xd(:)) ./ floor(max(im_out_size) / 2);
        
        for n = 1:length(new_period)
            tmp_im = double(stimulus.images(:,:,new_period(n)));
            tmp_im = imresize(tmp_im,im_out_size,'nearest');
            meg_stim.full_im(:,:,n) = double(ceil(abs(tmp_im - bk) ./ max(abs(stim_range - bk)))).*(srate.^2);
            
        end
        
        [meg_stim.full_x, meg_stim.full_y] = meshgrid(linspace(x_range(1), x_range(2), im_out_size(1)),...
            linspace(y_range(1), y_range(2), im_out_size(2)));
        
        keep_window = sum(meg_stim.full_im,3) > 0;
        
        meg_stim.im = reshape(meg_stim.full_im,[],size(meg_stim.full_im,3));
        meg_stim.im =  meg_stim.im(keep_window(:),:);
        meg_stim.X = meg_stim.full_x(:);
        meg_stim.Y = meg_stim.full_y(:);
        meg_stim.X = meg_stim.X(keep_window);
        meg_stim.Y = meg_stim.Y(keep_window);
        
        
        save(dest_file,'meg_stim');
        
        if exist(dest_file,'file')
            success = true;
            FULL_STIM_PATH = dest_file;
        end
        
        
        
        
    otherwise
        fprintf('Error: unknow stimulus type: %s\n',w_stim)
        
        
        
        
        
end


