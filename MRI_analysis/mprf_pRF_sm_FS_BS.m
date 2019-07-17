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

bs_model_path = dirPth.bs.dataPth;
bs_model_file = dir(fullfile(bs_model_path,'R*','headmodel_surf_os_meg.mat'));

bs_anat_path = dirPth.bs.anatPth;
bs_anat_file = strcat(bs_anat_path,'/subjectimage_T1.mat');


if ~exist(prf_dir_BS, 'dir')
    mkdir(prf_dir_BS);
    mkdir(roi_dir_BS);
end
% ----------

%% ------------------------------------------------------------------------
% Exporting pRF parameters from freesurfer vertices to brainstorm vertices
% -------------------------------------------------------------------------

load(fullfile(bs_model_file.folder,bs_model_file.name));
bs_head_model_surf = SurfaceFile;
% Load the anatomy imported by brainstorm. This file has a transformation
% maxtrix that is applied to the freesurfer surfaces when they are imported
% by brainstorm, so we need to apply the same transformation to the surface
% data as well:
bs_mri = load(bs_anat_file);

hs_to_load = {'lh','rh'};        % In this order as BS concatenates the hemispheres like this

% Check surface:
if any(strcmpi('pial',strsplit(bs_head_model_surf,'_')))
    surfaces_to_load = {'pial'};
    
    % (EK): Don't load white surface, since order of vertices is different
    % and will not match the downsampled prf data or rois
    % elseif any(strcmpi('white',strsplit(bs_head_model_surf,'_')))
    %     surfaces_to_load = {'white'};
    %
else
    error('Could not determine the surface to load');
    
end

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
    
    % Get brainstorm mesh
    surf_path = fullfile(bs_anat_path);
    
    % Load the brainstorm surface file:
    bs_surf_files = dir(fullfile(surf_path,'*tess_cortex_pial_low.mat'));
    
    load(fullfile(surf_path,bs_surf_files.name),'Vertices', 'Faces');
    
    bs_vert_idx = dsearchn(fs_vertices, Vertices);
    
    %% --------------------------------------------------------------------
    % Export prfs    
    % ---------------------------------------------------------------------
    
    pname = prf_dir_FS;
    
    % ROIs
    lh_files = dir(fullfile(pname,'lh.*'));
    
    % Loop over the LHs files (i.e. all the parameters on the lh surfaces):
    for n_hs = 1:length(lh_files)
        
        cur_lh_file = lh_files(n_hs).name;
        [~,~,par_name] = fileparts(cur_lh_file);
        
        par_name = par_name(2:end);
        % Find the corresponding rh file:
        cur_rh_file = ['rh.' par_name];
        
        % load and concatenate:
        both_data = [read_curv(fullfile(pname,cur_lh_file));...
            read_curv(fullfile(pname,cur_rh_file))];
        
        % preallocate the output variable:
        %             both_bs_data_out = nan(size(bs_vert_idx));
        %             both_bs_data_out = nan(size(bs_vert_idx));
        
        % Select to correct parameters:
        %             both_bs_data_out(bs_vert_idx) = both_data(bs_vert_idx);
        both_bs_data_out = both_data(bs_vert_idx);
        
        
        % Store the results:
        cur_out_file = [surfaces_to_load{n_surf} '.' par_name];
        
        fname = fullfile(prf_dir_BS,cur_out_file);
        
        write_curv(fname,both_bs_data_out,1);
        fprintf('%s\n',cur_out_file);
        
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
            
            cur_lh_rois = MRIread(fullfile(pname,cur_lh_file));
            cur_rh_rois = MRIread(fullfile(pname,cur_rh_file));
            
            both_data = [squeeze(cur_lh_rois.vol);squeeze(cur_rh_rois.vol)];
            
            both_bs_data_out = both_data(bs_vert_idx);
            both_bs_data_out_mask = both_bs_data_out;
            both_bs_data_out_mask(both_bs_data_out_mask > 0) = 1;
            
            % Store the results:
            cur_out_file = [surfaces_to_load{n_surf} '.' r];
            
            fname = fullfile(prf_dir_BS,cur_out_file);
            
            write_curv(fname,both_bs_data_out,1);
            fprintf('%s\n',cur_out_file);
            
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
            both_data = [read_curv(fullfile(pname,cur_lh_file));...
                read_curv(fullfile(pname,cur_rh_file))];
            
            % preallocate the output variable:
            %             both_bs_data_out = nan(size(bs_vert_idx));
            %             both_bs_data_out = nan(size(bs_vert_idx));
            
            % Select to correct parameters:
            %             both_bs_data_out(bs_vert_idx) = both_data(bs_vert_idx);
            both_bs_data_out = both_data(bs_vert_idx);
            
            % Store the results:
            cur_out_file = [surfaces_to_load{n_surf} '.' par_name];
            
            fname = fullfile(roi_dir_BS,cur_out_file);
            
            write_curv(fname,both_bs_data_out,1);
            fprintf('%s\n',cur_out_file);
        end
        
    end
    
    
end