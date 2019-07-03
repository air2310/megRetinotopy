function mprf_pRF_sm_FS_BS(dirPth)


% File paths
% ----------
rootDir = dirPth.rootPth;

freesurfer_surface = strcat(rootDir,dirPth.freesurfer.surfPth(2:end));
prf_dir_FS = strcat(rootDir,dirPth.fmri.saveDataPth_prfFS(2:end));
roi_dir_FS = strcat(rootDir,dirPth.fmri.saveDataPth_roiFS(2:end));

prf_dir_BS = strcat(rootDir,dirPth.fmri.saveDataPth_prfBS(2:end));
roi_dir_BS = strcat(rootDir,dirPth.fmri.saveDataPth_roiBS(2:end));


bs_model_path = strcat(rootDir,dirPth.brainstorm.dataPth(2:end));
bs_model_file = strcat(bs_model_path,'/headmodel_surf_os_meg.mat');
        
bs_anat_path = strcat(rootDir,dirPth.brainstorm.anatPth(2:end));
bs_anat_file = strcat(bs_anat_path,'/subjectimage_T1.mat');

% ----------

load(bs_model_file);
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

for n=1:length(surfaces_to_load)

    % Get brainstorm mesh
    fs_vertices = [];
    % Loop over all (2, left and right) surfaces:
    for nn = 1:length(hs_to_load)
        cur_surf = [hs_to_load{nn} '.' surfaces_to_load{n}];
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
        verts = [fs_vertices; tmp_verts];  
    end
     
    % Get brainstorm mesh
    surf_path = fullfile(bs_anat_path);
   
    % Load the brainstorm surface file:
    bs_surf_files = dir(fullfile(surf_path,'*tess_cortex_pial_low.mat'));
    
    load(fullfile(surf_path,bs_surf_files.name),'Vertices', 'Faces');
        
    bs_vert_idx = dsearchn(verts, Vertices);
    
    %%%%%%%%%%%%%%%%%%
    % Export prfs    %
    %%%%%%%%%%%%%%%%%%
    
    pname = prf_dir_FS;
    w_export = 'prf';
    
    % ROIs 
    lh_files = dir(fullfile(pname,'lh.*'));    
        
    % Loop over the LHs files (i.e. all the parameters on the lh surfaces):
    for nn = 1:length(lh_files)
        
        cur_lh_file = lh_files(nn).name;
        [~,~,par_name] = fileparts(cur_lh_file);
        if strcmpi(par_name,'.mgz')
            
            
        else
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
            cur_out_file = [surfaces_to_load{n} '.' par_name];
            
            if strcmpi(w_export,'prf')
                fname = fullfile(prf_dir_BS,cur_out_file);
                
            elseif strcmpi(w_export,'roi')
                fname = fullfile(roi_dir_BS,cur_out_file);
                
            end
            
            write_curv(fname,both_bs_data_out,1);
            fprintf('%s\n',cur_out_file);
        end
        
    end    
    
    
    
    %%%%%%%%%%%%%%%%%%
    % Export rois    %
    %%%%%%%%%%%%%%%%%%
    
    pname = roi_dir_FS;
    w_export = 'roi';
    
    % ROIs 
    lh_files = dir(fullfile(pname,'lh.*'));    
        
    % Loop over the LHs files (i.e. all the parameters on the lh surfaces):
    for nn = 1:length(lh_files)
        
        cur_lh_file = lh_files(nn).name;
        [~,~,par_name] = fileparts(cur_lh_file);
        if strcmpi(par_name,'.mgz')
            
            
        else
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
            cur_out_file = [surfaces_to_load{n} '.' par_name];
            
            if strcmpi(w_export,'prf')
                fname = fullfile(prf_dir_BS,cur_out_file);
                
            elseif strcmpi(w_export,'roi')
                fname = fullfile(roi_dir_BS,cur_out_file);
                
            end
            
            write_curv(fname,both_bs_data_out,1);
            fprintf('%s\n',cur_out_file);
        end
        
    end    

        
end