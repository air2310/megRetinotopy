function mprf_pRF_sm_FS_BS(dirPth,opt)
% mprf_pRF_sm_FS_BS(dirPth,plot_stim)
%
% Function to export prf parameters (unsmooth/smooth) from mrVista Gray
% Ribbon (volume) to freesurfer surface vertices
%                  
% INPUTS:
%   dirPth      :   paths locating subject's data and files (struct, see loadPaths.m)
%   opt         :   struct with boolean flags, needed to request using rois
%                   defined in mrVista space 

%% ---------- 
% File paths
% -----------
freesurfer_surface = dirPth.fs.surfPth;
prf_dir_FS = dirPth.fmri.saveDataPth_prfFS;
roi_dir_FS = dirPth.fmri.saveDataPth_roiFS;
prf_dir_BS = dirPth.fmri.saveDataPth_prfBS;
roi_dir_BS = dirPth.fmri.saveDataPth_roiBS;

bs_anat_path = dirPth.bs.anatPth;
bs_T1anat_file = strcat(bs_anat_path,'/subjectimage_T1.mat');

if ~exist(prf_dir_BS, 'dir')
    mkdir(prf_dir_BS);
    mkdir(roi_dir_BS);
end
% ----------

%% ------------------------------------------------------------------------
% Exporting pRF parameters from freesurfer vertices to brainstorm vertices
% -------------------------------------------------------------------------

% Load the anatomy imported by brainstorm. This file has a transformation
% maxtrix that is applied to the freesurfer surfaces when they are imported
% by brainstorm, so we need to apply the same transformation to the surface
% data as well:
bs_mri = load(bs_T1anat_file);

hs_to_load       = {'lh','rh'};        % In this order as BS concatenates the hemispheres like this
surfaces_to_load = {'pial'};


fprintf('Exporting surface data:\n');

for n_surf=1:length(surfaces_to_load)
    
    % Get brainstorm mesh
    fs_vertices = [];
    % Loop over all (2, left and right) surfaces:
    for n_hs = 1:length(hs_to_load)
        cur_surf = [hs_to_load{n_hs} '.' surfaces_to_load{n_surf}];
        surf_file = fullfile(freesurfer_surface,cur_surf);
        
        % Use the same routine as Brainstorm uses, otherwise the vertices
        % are slightly off and intersectCols does not find all the matching
        % vertices.
        
        [tmp_verts, ~] = mne_read_surface(surf_file);
        
        % Transformations applied by brainstorm:
        tmp_verts = bsxfun(@plus, tmp_verts, [128 129 128] / 1000);
        tmp_verts = [tmp_verts'; ones(1,size(tmp_verts,1))]; %% Actually verts)
        tmp_verts = [bs_mri.SCS.R, bs_mri.SCS.T./1000; 0 0 0 1] * tmp_verts;
        tmp_verts = tmp_verts(1:3,:)';
        
        
        % combine left and right hemi
        fs_vertices = [fs_vertices; tmp_verts];
    end
    
    % Load the brainstorm pial surface file:
    bs_pial_surf_file = dir(fullfile(bs_anat_path,'*tess_cortex_pial_low.mat'));
    
    load(fullfile(bs_anat_path,bs_pial_surf_file.name),'Vertices', 'Faces');
    
    bs_vert_idx = dsearchn(fs_vertices, Vertices);
    
    %% --------------------------------------------------------------------
    % Export prfs    
    % ---------------------------------------------------------------------
    
    pname = prf_dir_FS;
    
    lh_files = dir(fullfile(pname,'lh.*'));
    
    for nn = 1:length(lh_files)
    
        cur_lh_file = lh_files(nn).name;
        [~,~,par_name] = fileparts(cur_lh_file);
        par_name = strsplit(par_name, '.');
        par_name = par_name{2};
    
        prfParamfname = @(hem)(sprintf('%s/%s.%s', pname, hem, par_name));
    
        surfdataFS.lh = read_curv(prfParamfname('lh')); 
        surfdataFS.rh = read_curv(prfParamfname('rh'));
    
    
        % load and concatenate:
        both_fs_data = [surfdataFS.lh; surfdataFS.rh];
        
        % Find corresponding BS vertices for FS surface data
        both_bs_data_out = both_fs_data(bs_vert_idx);
        
        % Store the BS results:
        cur_out_file = [surfaces_to_load{n} '.' par_name]; 
        fname = fullfile(prf_dir_BS,cur_out_file);
        write_curv(fname,both_bs_data_out,1);
        fprintf('(%s): Brainstorm combined hemi files: %s\n',mfilename, cur_out_file);
        
        % Store the FS results:
        cur_out_file = [surfaces_to_load{n} '.' par_name]; 
        fname = fullfile(prf_dir_FS,cur_out_file);
        write_curv(fname,both_fs_data,1);
        fprintf('(%s): Freesurfer combined hemi files: %s\n',mfilename, cur_out_file);

        
    end
    
    
    %% --------------------------------------------------------------------
    % Export rois
    %----------------------------------------------------------------------
    if opt.roimrvToFS == 1
        pname = roi_dir_FS;
        lh_files = dir(fullfile(pname,'lh.*')); % ROIs
    else
        pname = dirPth.fs.surfPth;
        lh_files = dir(fullfile(pname,'lh.wang2015_atlas.mgz')); % ROIs
    end
    
    % Loop over the LHS files (i.e. all the parameters on the lh surfaces):
    for n_hs = 1:length(lh_files)
        
        cur_lh_file = lh_files(n_hs).name;
        [~,r,par_name] = fileparts(cur_lh_file);
        
        if strcmpi(par_name,'.mgz') % To load Wang et al rois
            
            r = r(4:end);
            
            % Find the corresponding rh file:
            cur_rh_file = ['rh.',r, par_name];
            
            % Read in both hemispheres' vertices
            cur_lh_rois = MRIread(fullfile(pname,cur_lh_file));
            cur_rh_rois = MRIread(fullfile(pname,cur_rh_file));
            
            % Combine the hemispheres
            both_fs_data = [squeeze(cur_lh_rois.vol);squeeze(cur_rh_rois.vol)];
            
            % Resample to nearest Brainstorm vertices
            both_bs_data_out = both_fs_data(bs_vert_idx);
            
            
            % Store full Wang atlas as BS file
            cur_out_file = [surfaces_to_load{n} '.' r];           
            fname = fullfile(roi_dir_BS,cur_out_file); 
            write_curv(fname,both_bs_data_out,1);
            fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, cur_out_file);
            
            % Store combined Wang atlas as FS file           
            fname = fullfile(roi_dir_FS,cur_out_file); 
            write_curv(fname,both_fs_data,1);
            fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, cur_out_file);
            
            % Create roi mask BS pial file
            BS_mask = both_bs_data_out>0;
            cur_out_file = [surfaces_to_load{n} '.mask'];
            fname = fullfile(prf_dir_BS,cur_out_file);
            write_curv(fname,BS_mask,1);
            fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, cur_out_file);
            
            % Create roi mask FS pial file
            FS_mask = both_bs_data_out>0;
            fname = fullfile(prf_dir_FS,cur_out_file);
            write_curv(fname,FS_mask,1);
            fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, cur_out_file);

            
            % Create V1-V3 Wang roi BS mask pial file
            V123idx = 1:6; % first six rois are "V1v" "V1d" "V2v" "V2d" "V3v" "V3d"
            
            BS_V123mask = ismember(both_bs_data_out, V123idx);
            cur_out_file = [surfaces_to_load{n} '.V123mask'];
            fname = fullfile(prf_dir_BS,cur_out_file);
            write_curv(fname,BS_V123mask,1);
            fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, cur_out_file);
            
            % Create V1-V3 Wang roi FS mask lh, rh file
            FS_V123mask = ismember(both_fs_data, V123idx);
            fname = fullfile(prf_dir_FS,cur_out_file);
            write_curv(fname,FS_V123mask,1);
            fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, cur_out_file);

            
                
            
            % store the result for the mask
            cur_out_file_mask = [surfaces_to_load{n_surf},'.','mask'];
            
            fname_mask = fullfile(prf_dir_BS,cur_out_file_mask); % save the mask in the pRF directory instead of the roi directory
            
            write_curv(fname_mask,both_bs_data_out_mask,1);
            fprintf('%s\n',cur_out_file_mask);
            
            
        else % Load ROIs drawn in mrVista surface and exported to freesurfer space
            par_name = par_name(2:end);
            % Find the corresponding rh file:
            cur_rh_file = ['rh.' par_name];
            
            % load and concatenate:
            both_fs_data = [read_curv(fullfile(pname,cur_lh_file));...
                read_curv(fullfile(pname,cur_rh_file))];
            
            % preallocate the output variable:
            %             both_bs_data_out = nan(size(bs_vert_idx));
            %             both_bs_data_out = nan(size(bs_vert_idx));
            
            % Select to correct parameters:
            %             both_bs_data_out(bs_vert_idx) = both_data(bs_vert_idx);
            both_bs_data_out = both_fs_data(bs_vert_idx);
            
            % Store the results:
            cur_out_file = [surfaces_to_load{n} '.' par_name];        
            fname = fullfile(roi_dir_BS,cur_out_file);  
            write_curv(fname,both_bs_data_out,1);
            fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, cur_out_file);
            
            
            % Create roi mask FS pial file
            fname = fullfile(prf_dir_FS,cur_out_file);
            write_curv(fname,both_fs_data,1);
            fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, cur_out_file);

            
            
        end
        
    end
    
    
end