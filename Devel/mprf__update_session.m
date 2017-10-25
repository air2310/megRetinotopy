function mprf__update_session(do_stim)
% Called by:
%   - mprfSession_init
% Check if we want to do the stimulus here, or later on. I do not know why
% this is here, but the stimulus is imported later on by the
% mprfSession_init routine.
if ~exist('do_stim','var') || isempty(do_stim)
    do_stim.rm_stim = false;
    do_stim.meg_stim = false;
end

% Declare mprfSESSION as global. It Assumes that a 
global mprfSESSION

% The mprf__session_init_gui routine adds an 'updated' field to the
% mprfSESSION variable. This field marks which paths have changed and so
% which part of the session need updating. 
if isfield(mprfSESSION,'updated')
    % Get all the fields from the updated field:
    to_update = getAllFieldsFromStruct(mprfSESSION.updated);
    
else
    fprintf('Nothing to update\n')
    return
end

% Loop over the fields that need updating:
for n = 1:length(to_update)
    fname = '';
    fext = '';
    cur_update = to_update{n};
    % If an original paths needs updating:
    if any(strfind(cur_update,'.orig.paths'))
        orig_path = eval(['mprfSESSION' cur_update]);
        
        % Update original path and import data again. Generally, this
        % follows the following simple logic:
        % 1. Create the source path (orig_path)
        % 2. Determine the file type
        % 3. create the output file name
        % 4. set the yes to all variable to false
        
        % This becomes the input for the mprf__import_data function. This
        % function checks if the already exists and asks if you want to
        % overwrite the existing data (if yes_to_all = true, then it wont
        % bother you with this, it just overwrites it).
        % mprf__import_data copies the data and returns the yes_to_all
        % variable, as changes when you request it through the dialogue. It
        % also returns a success variable, which tells you if the output
        % filename is indeed created or not.
        % For some data types it set the file name in the mprfSESSION
        % variable by using mprf__session_set_file. This way,
        % mprf__session_get_file can be used to get the file directly.
        
        if (any(strfind(cur_update, 'vista_t1')) || any(strfind(cur_update, 'vista_class')))
            
            if strfind(cur_update,'vista_t1')
                cur_vol = 'vista_t1';
            elseif strfind(cur_update,'vista_class')
                cur_vol = 'vista_class';
            end
            
            yes_to_all = false;
            fname = orig_path(find(orig_path == '/',1,'last')+1 : find(orig_path == '.',1,'first')-1);
            fext = orig_path(find(orig_path == '.',1,'first') : end);
            
            
            [yes_to_all,succes] = mprf__import_data(orig_path,cur_vol,[fname, fext],yes_to_all);
            
            if succes
                mprfSESSION.has.(cur_vol) = true;
            else
                mprfSESSION.has.(cur_vol) = false;
                
            end

            mprf__set_file(cur_vol, [fname, fext]);
            
        elseif (any(strfind(cur_update, 'lh_rois')) || ...
                any(strfind(cur_update, 'rh_rois')))
            
            if strfind(cur_update, 'lh_rois')
                cur_rois = 'lh_rois';
            elseif strfind(cur_update, 'rh_rois')
                cur_rois = 'rh_rois';
            end
            
            if iscell(orig_path)
            else
                orig_path = {orig_path};
                
            end
            
            succes = zeros(size(orig_path));
            yes_to_all = false;
            for nn = 1:length(orig_path)
                cur_orig_path = orig_path{nn};
                
                [~,fname,fext] = fileparts(cur_orig_path);
                
                [yes_to_all, succes(nn)] = mprf__import_data(cur_orig_path, cur_rois,[fname,fext] ,yes_to_all);

            end
            
            if any(succes)
                mprfSESSION.has.(cur_rois) = true;
                
            else
                mprfSESSION.has.(cur_rois) = false;
                
            end
                        
        elseif any(strfind(cur_update, 'meg_stim_params'))
            indir = dir(fullfile(orig_path,'*.mat'));
            
            succes = zeros(size(indir));
            yes_to_all = false;
            
            for nn = 1:length(indir)
                cur_orig_path = fullfile(orig_path,indir(nn).name);
                
                [~,fname,fext] = fileparts(cur_orig_path);
                
                [yes_to_all, succes(nn)] = mprf__import_data(cur_orig_path, 'meg_stim_params',[fname,fext] ,yes_to_all);
                
            end
            
            if any(succes)
                mprfSESSION.has.meg_stim_params = true;
                
            else
                mprfSESSION.has.meg_stim_params = false;
                
            end
            
            
            
        elseif any(strfind(cur_update, 'meg_stim'))
            
            indir = dir(fullfile(orig_path,'*.mat'));
            
            succes = zeros(size(indir));
            yes_to_all = false;
            
            for nn = 1:length(indir)
                cur_orig_path = fullfile(orig_path,indir(nn).name);
                
                [~,fname,fext] = fileparts(cur_orig_path);
                
                [yes_to_all, succes(nn)] = mprf__import_data(cur_orig_path, 'meg_stim',[fname,fext] ,yes_to_all);
                
            end
            
            if any(succes)
                mprfSESSION.has.meg_stim = true;
                
            else
                mprfSESSION.has.meg_stim = false;
                
            end
            
            
            
        elseif any(strfind(cur_update, 'mri_stim'))
            indir = dir(fullfile(orig_path,'*.mat'));
            
            succes = zeros(size(indir));
            yes_to_all = false;
            
            for nn = 1:length(indir)
                cur_orig_path = fullfile(orig_path,indir(nn).name);
                
                [~,fname,fext] = fileparts(cur_orig_path);
                
                [yes_to_all, succes(nn)] = mprf__import_data(cur_orig_path, 'mri_stim',[fname,fext] ,yes_to_all);
                
            end
            
            if any(succes)
                mprfSESSION.has.mri_stim = true;
                
            else
                mprfSESSION.has.mri_stim = false;
                
            end
            
            
            
            
        elseif any(strfind(cur_update, 'meg_data'))
            yes_to_all = false;
            [~,fname,fext] = fileparts(orig_path);
            
            [yes_to_all, succes] = mprf__import_data(orig_path, 'meg_data',[fname,fext] ,yes_to_all);
            
            if succes
                mprfSESSION.has.meg_data = true;
                
            else
                mprfSESSION.has.meg_data = false;
                
            end
            
            
            
            
        elseif any(strfind(cur_update, 'rm_model'))
            yes_to_all = false;
            [~,fname,fext] = fileparts(orig_path);
            
            [yes_to_all, succes] = mprf__import_data(orig_path, 'rm_model',[fname,fext] ,yes_to_all);
            
            if succes
            mprfSESSION.has.rm_model = true;
            
            else
                mprfSESSION.has.rm_model = false;
                
            end
            
            
        elseif any(strfind(cur_update, 'bs_hm'))
            yes_to_all = false;
            load(orig_path,'Gain','GridOrient','GridLoc','SurfaceFile')
            
            [~,this_surf] = fileparts(SurfaceFile);
            out_name = ['head_model_' this_surf '.mat'];
            dest_path = fullfile(mprf__get_directory('main_dir'),mprf__get_directory('bs_headmodel'),out_name);
            
            bs_model.Gain = Gain;
            bs_model.GridOrient = GridOrient;
            bs_model.SurfaceFile = SurfaceFile;
            bs_model.GridLoc = GridLoc;
            
            if exist(dest_path,'file')
                
                answer = questdlg(sprintf('%s is already imported, overwrite?',out_name),...
                    'Overwrite','No','Yes','Yes to all','No');
                
                switch lower(answer)
                    
                    case 'no'
                        
                    case {'yes', 'yes to all'}
                        save(dest_path,'bs_model');
                        
                end
            else
                save(dest_path,'bs_model');
                
            end
            
            if exist(dest_path,'file')
                fprintf('Imported brainstorm head model file\n')
                mprfSESSION.has.bs_hm = true;
                
            else
                mprfSESSION.has.bs_hm = false;
                
                
            end
            
            
            
        elseif any(strfind(cur_update, 'fs_dir'))
            succes = 0;
            if exist(fullfile(orig_path,'mri'),'dir');
                if exist(fullfile(orig_path,'mri','t1.mgz'),'file');
                    yes_to_all = false;
                    cur_orig_path = fullfile(orig_path,'mri','t1.mgz');
                    
                    [~,rm_fname,rm_ext] = fileparts(orig_path);
                    
                    [yes_to_all, succes] = mprf__import_data(cur_orig_path, 'fs_volume','t1.mgz' ,yes_to_all);

                else
                    fprintf('Could not locate freesurfer T1 in mri directory\n')
                end
                
            else
                fprintf('Could not locate freesurfer mri directory: %s \n',fullfile(orig_path,'mri'))
                
            end
            
            if succes
                mprfSESSION.has.fs_volume = true;
                
            else
                mprfSESSION.has.fs_volume = false;
                
            end
            
            surfaces_to_import = {'lh.pial','rh.pial','lh.white','rh.white'};
            
            if exist(fullfile(orig_path,'surf'),'dir');
                yes_to_all = false;
                
                
                for nn = 1:length(surfaces_to_import)
                    cur_surf = surfaces_to_import{nn};
                    
                    cur_orig_path = fullfile(orig_path, 'surf',cur_surf);
                    succes = 0;
                    if exist(cur_orig_path,'file')
                        [yes_to_all, succes] = mprf__import_data(cur_orig_path,'fs_surface',cur_surf,yes_to_all);
                        
                    else
                        
                    end
                    cur_surf(cur_surf == '.') = '_';
                    [~, tmp_name] = fileparts(cur_surf);
                    if succes
                        mprfSESSION.has.fs_surface.(tmp_name) = true;
                        
                    else
                        mprfSESSION.has.fs_surface.(tmp_name) = false;
                        
                    end
                    
                    
                end
                
            else
                
            end

        elseif any(strfind(cur_update, 'bs_anat'))
            
            
            files_to_import = {'subjectimage_T1.mat','tess_cortex_mid_low.mat',...
                'tess_cortex_pial_low.mat','tess_cortex_white_low.mat'};
            
            yes_to_all = false;
            
            for nn = 1:length(files_to_import)
                cur_file = files_to_import{nn};
                cur_orig_path = fullfile(orig_path, cur_file);
                succes = 0;
                if exist(fullfile(orig_path, cur_file),'file');
                    [yes_to_all, succes] = mprf__import_data(cur_orig_path,'bs_anat',cur_file,yes_to_all);
                end
                [~, tmp_name] = fileparts(cur_file);
                
                if succes
                    mprfSESSION.has.bs_anat.(tmp_name) = true;
                    
                else
                    mprfSESSION.has.bs_anat.(tmp_name) = false;
                    
                end
                
                
            end

            mprf__set_file('bs_anat','subjectimage_T1.mat')
            
            
        elseif any(strfind(cur_update, 'vista_path'));
            
            if exist(mprfSESSION.orig.paths.vista_path,'dir')
                mprfSESSION.has.vista_path = true;
                
            else
                mprfSESSION.has.vista_path = false;
                
            end
   
        end
        
    end

    if n == length(to_update)
        mprfSESSION = rmfield(mprfSESSION,'updated');
        
        
    end
end




end








