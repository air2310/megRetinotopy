function mprfSessionImportSourceData(fields)
global mprfSESSION

if ~exist(fullfile(mprfSESSION.init.main_dir,'source'),'dir')
    mkdir(pwd,'source')
end

if ~exist('fields','var') || isempty(fields) % Assume to loop over all orig fields in the mprfSESSION struct
    fields = fieldnames(mprfSESSION.orig);
    loop_over_paths = true;
    incr = 1;
    
else
    loop_over_paths = false;
    incr = 2;
end

for n = 1:incr:length(fields)
    
    
    cur_field = fields{n};
    
    if loop_over_paths
        orig_path = mprfSESSION.orig.(cur_field);
    end
    
    if strcmpi(cur_field,'ret_model') % Also import T1 and class file from mrVista
        
        if ~exist(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source'),'mrvista')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'ret_model')
        elseif ~exist(fullfile(mprfSESSION.init.main_dir,'source','mrvista','ret_model'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'ret_model')
        end
        
        [~,rm_fname,rm_ext] = fileparts(mprfSESSION.orig.(cur_field));
        
        if exist(fullfile(mprfSESSION.init.main_dir,'source','mrvista','ret_model',[rm_fname rm_ext]),'file')
            fprintf('Found retinotopy file with the same name, skipping import\n')
            continue
        end
        dest_path = fullfile(mprfSESSION.init.main_dir,'source','mrvista','ret_model',[rm_fname rm_ext]);
        
        sys_cmnd = sprintf('cp %s %s',orig_path, dest_path);
        
        if system(sys_cmnd)
            warning('Failed to copy retinotopy model file\n')
        else
            mprfUpdateSession(3,'add_source_path','ret_model',dest_path);
            fprintf('Imported retinotopy model file\n')
        end
        
        
    elseif strcmpi(cur_field,'bs_head_model')
        
        if ~exist(fullfile(mprfSESSION.init.main_dir,'source','brainstorm'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source'),'brainstorm')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','brainstorm'),'head_model')
        elseif ~exist(fullfile(mprfSESSION.init.main_dir,'source','brainstorm','head_model'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','brainstorm'),'head_model')
        end
        
        
        load(orig_path,'Gain','GridOrient','GridLoc','SurfaceFile')
        
        [~,this_surf] = fileparts(SurfaceFile);
        out_name = ['head_model_' this_surf '.mat'];
        dest_path = fullfile(mprfSESSION.init.main_dir,'source','brainstorm','head_model',out_name);
        
        if exist(dest_path,'file')
            fprintf('Found Brainstorm model file with the same name, skipping import\n')
            continue
        end
        
        bs_model.Gain = Gain;
        bs_model.GridOrient = GridOrient;
        bs_model.SurfaceFile = SurfaceFile;
        bs_model.GridLoc = GridLoc;
        
        save(dest_path,'bs_model');
        mprfUpdateSession(3,'add_source_path','bs_head_model',dest_path);
        fprintf('Imported brainstorm head model file\n')
        
        
    elseif strcmpi(cur_field,'fs_subject_dir')
        
        if exist(fullfile(orig_path,'mri'),'dir');
            if exist(fullfile(orig_path,'mri','t1.mgz'),'file');
                
                if ~exist(fullfile(mprfSESSION.init.main_dir,'source','freesurfer'),'dir')
                    mkdir(fullfile(mprfSESSION.init.main_dir,'source'),'freesurfer')
                    mkdir(fullfile(mprfSESSION.init.main_dir,'source','freesurfer'),'volumes')
                elseif ~exist(fullfile(mprfSESSION.init.main_dir,'source','freesurfer','volumes'),'dir')
                    mkdir(fullfile(mprfSESSION.init.main_dir,'source','freesurfer'),'volumes')
                end
                
                dest_path = fullfile(mprfSESSION.init.main_dir,'source','freesurfer','volumes','t1.mgz');
                
                if exist(dest_path,'file')
                    fprintf('Freesurfer T1 already imported, skipping\n')
                    
                else
                    
                    sys_cmnd = sprintf('cp %s %s',fullfile(orig_path,'mri','t1.mgz'),...
                        dest_path);
                    
                    if system(sys_cmnd)
                        warning('Failed to copy freesurfer T1 file\n')
                    else
                        mprfUpdateSession(3,'add_source_path','fs_t1',dest_path);
                        fprintf('Imported freesurfer t1\n')
                        
                    end
                end
            else
                fprintf('Could not locate freesurfer T1 in mri directory\n')
            end
            
            
        else
            fprintf('Could not locate freesurfer mri directory: %s \n',fullfile(orig_path,'mri'))
            
            
        end
        
        if exist(fullfile(orig_path,'surf'),'dir');
            
            if ~exist(fullfile(mprfSESSION.init.main_dir,'source','freesurfer'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'source'),'freesurfer')
                mkdir(fullfile(mprfSESSION.init.main_dir,'source','freesurfer'),'surfaces')
            elseif ~exist(fullfile(mprfSESSION.init.main_dir,'source','freesurfer','surfaces'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'source','freesurfer'),'surfaces')
            end
            
            surfaces_to_import = {'lh.pial','rh.pial','lh.white','rh.white'};
            
            for nn = 1:length(surfaces_to_import)
                cur_surf = surfaces_to_import{nn};
                dest_path = fullfile(mprfSESSION.init.main_dir,'source','freesurfer','surfaces',cur_surf);
                if exist(fullfile(orig_path,'surf',cur_surf),'file') % exist(fullfile(orig_path,'surf','t1.mgz'),'file'); %
                    
                    if exist(dest_path,'file')
                        fprintf('Freesurfer surface %s already imported, skipping\n',cur_surf)
                        continue
                        
                    else
                        
                        sys_cmnd = sprintf('cp %s %s',fullfile(orig_path,'surf',cur_surf),...
                            dest_path);
                        
                        if system(sys_cmnd)
                            warning('Failed to copy freesurfer %s surface\n',cur_surf)
                        else
                            cur_field_name = cur_surf;
                            cur_field_name(cur_field_name == '.') = '_';
                            mprfUpdateSession(3,'add_source_path',cur_field_name,dest_path);
                            fprintf('Imported freesurfer %s surface\n',cur_surf)
                            
                            
                        end
                        
                    end
                else
                    fprintf('Could not find freesurfer %s surface\n',cur_surf)
                    
                end
                
            end
            mprfUpdateSession(3,'add_source_path','fs_surface_dir',fileparts(dest_path));
            
            
        else
            fprintf('Could not locate freesurfer surface directory: %s \n',fullfile(orig_path,'surf'))
            
            
        end
        
        
        
    elseif strcmpi(cur_field,'bs_anat_dir')
        
        if ~exist(fullfile(mprfSESSION.init.main_dir,'source','brainstorm'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source'),'brainstorm')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','brainstorm'),'anatomy')
        elseif ~exist(fullfile(mprfSESSION.init.main_dir,'source','brainstorm','anatomy'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','brainstorm'),'anatomy')
        end
        
        
        files_to_import = {'subjectimage_T1.mat','tess_cortex_mid_low.mat',...
            'tess_cortex_pial_low.mat','tess_cortex_white_low.mat'};
        
        
        
        for nn = 1:length(files_to_import)
            
            cur_file = files_to_import{nn};
            dest_path = fullfile(mprfSESSION.init.main_dir,'source','brainstorm','anatomy',cur_file);
            if exist(fullfile(orig_path,cur_file),'file') % exist(fullfile(orig_path,'surf','t1.mgz'),'file'); %
                
                if exist(dest_path,'file')
                    fprintf('Brainstorm file %s already imported, skipping\n',cur_file)
                    continue
                    
                end
                
                sys_cmnd = sprintf('cp %s %s',fullfile(orig_path,cur_file),...
                    dest_path);
                
                if system(sys_cmnd)
                    warning('Failed to copy brainstorm %s file',cur_file)
                else
                    [~,cur_field_name] = fileparts(cur_file);
                    mprfUpdateSession(3,'add_source_path',cur_field_name,dest_path);
                    fprintf('Imported brainstorm %s file\n',cur_file)
                    
                    
                end
                
                
            else
                fprintf('Could not locate %s in Brainstorm anatomy directory',cur_file)
            end
        end
        
        
    elseif strcmpi(cur_field,'mr_vista_t1')
        
        if ~exist(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source'),'mrvista')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'anatomy')
        elseif ~exist(fullfile(mprfSESSION.init.main_dir,'source','mrvista','anatomy'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'anatomy')
        end
        
        fname = orig_path(find(orig_path == '/',1,'last')+1 : find(orig_path == '.',1,'first')-1);
        f_ext = orig_path(find(orig_path == '.',1,'first') : end);
        
        if ischar(orig_path)
            dest_path = fullfile(mprfSESSION.init.main_dir,'source','mrvista','anatomy',[fname,f_ext]);
            
            sys_cmnd = sprintf('cp %s %s',orig_path,dest_path);
            if system(sys_cmnd);
                warning('Could not copy mrVista T1\n');
            else
                mprfUpdateSession(3,'add_source_path','mrv_T1',dest_path);
                fprintf('Imported mrVista T1\n')
                
            end
            
            
        else
            warning('No valid file name or path returned, skipping')
        end
        
        
        
    elseif strcmpi(cur_field,'mr_vista_class')
        
        if ~exist(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source'),'mrvista')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'anatomy')
        elseif ~exist(fullfile(mprfSESSION.init.main_dir,'source','mrvista','anatomy'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'anatomy')
        end
        
        
        fname = orig_path(find(orig_path == '/',1,'last')+1 : find(orig_path == '.',1,'first')-1);
        f_ext = orig_path(find(orig_path == '.',1,'first') : end);
        
        if ischar(orig_path)
            dest_path = fullfile(mprfSESSION.init.main_dir,'source','mrvista','anatomy',[fname,f_ext]);
            
            sys_cmnd = sprintf('cp %s %s',orig_path,dest_path);
            if system(sys_cmnd);
                warning('Could not copy mrVista CLASS\n');
            else
                mprfUpdateSession(3,'add_source_path','mrv_class',dest_path);
                fprintf('Imported mrVista CLASS\n')
                
            end
            
            
        else
            warning('No valid file name or path returned, skipping')
        end
        
        
    elseif strcmpi(cur_field,'vista_dir')
        
        % Do nothing
        
        
    elseif strcmpi(cur_field,'fssubjectsdir')
        
        % Do nothing
        
    elseif strcmpi(cur_field,'rm_stimulus')
        if ~exist(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source'),'mrvista')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'rm_stimulus')
        elseif ~exist(fullfile(mprfSESSION.init.main_dir,'source','mrvista','rm_stimulus'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'source','mrvista'),'rm_stimulus')
        end
        
        dest_path = fullfile(mprfSESSION.init.main_dir,'source','mrvista','rm_stimulus','rm_stim');
        
        if isstruct(fields{n+1})
            rm_stim = fields{n+1};
            save(dest_path,'rm_stim')
            fprintf('RM stimulus save as: %s\n',dest_path)
            mprfUpdateSession(3,'add_source_path','rm_stim',dest_path);
            
        else
           % Do not know yet... 
            
        end
    
    else
        warning('Unknown field %s. Skipping',cur_field)
        
        
        
    end
    
    
    
end


end