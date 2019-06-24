% Quality checks to see everything is working correctly
%
% 1) Projection of pRF parameters and ROIs from each voxel (mrVista space)
%    to each vertex (freesurfer space)
% 2) ROI projection / Predicted responses from each vertex (stimulus * pRF models)
% 3) Projection of parameters / ROIs onto brainstorm space
% 4) Synthetic dataset to see if pipeline works (Measured signal =
%    Predicted signal + noise ?)
% 5) Plot all 19 reference phases for a sensor with high VE
% 6) VE maps and predicted responses for every size and position range


% Just a hack for now, there should be a correct way of checking if the
% files are added to the path or not
if isempty(which('mprf_computemetric.m'))
    tbUse({'retMeg','vistasoft'});
end

%% for sub 004 - pRF parameters in mrVista space [Smooth the parameters (x,y,sigma, beta)]
sub = 'wlsubj068';
data_type = 'Averages';
rm_model = sprintf('/mnt/storage_2/MEG/Retinotopy/Data/fMRI/%s/vistaSession/Gray/Averages/rm_Averages-fFit.mat',sub);
seg_file = sprintf('/mnt/storage_2/MEG/Retinotopy/Data/Anatomy/%s/t1_class.nii.gz',sub);
mrVNif_prf_data = sprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/prf_data/nifti',sub);
rm_stim_file = sprintf('/mnt/storage_2/MEG/Retinotopy/Subject_sessions/%s/source/mrvista/rm_stimulus/rm_stim_rm_Averages-fFit.mat',sub);
prf_mat = sprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/prf_data/data_dir',sub);

% We need a volume view:

hvol = initHiddenGray;
anat_file = sprintf('/mnt/storage_2/MEG/Retinotopy/Data/Anatomy/%s/t1.nii.gz',sub);
setVAnatomyPath(anat_file);
anat_value = niftiRead(anat_file);
hvol = viewSet(hvol,'anatomy',anat_value.data);

% Set the volume view to the current data type and add the RM model
hvol = viewSet(hvol,'curdt',data_type);
hvol = rmSelect(hvol,1,rm_model);
% Mask to exclude unreliable voxels (i.e. VE == 0) from the smoothing
% below. Otherwise, the pRF parameters will be averaged with a lot of
% zeros:
sm_mask = rmGet(hvol.rm.retinotopyModels{1},'varexplained') > 0;
% We need these parameters from the pRF model
params = {'sigma','x','y','varexplained','beta'};
% We need the mrVista segmentation to check if the selection of pRF
% parameters is correct, i.e. all selected pRF parameters must fall in the
% gray matter.
cls = niftiRead(seg_file);

% stimulus file 
load(rm_stim_file);
%rm_stim =  rm_stim;

%%

% Initialize the weighted connectivity matrix used in the smoothing
wConMat = [];


for nn = 1:length(params)
    
    cur_param = params{nn};
    % current parameter's NIFTI file:
    fname = fullfile(mrVNif_prf_data,[cur_param '.nii.gz']);
    
    % Load the current parameter in the VOLUME view. This is mainly used to
    % export the pRF parameters to nifti files. The nifti files are really
    % just for inspection and are not used in any further analyses
    hvol = viewSet(hvol,'curdt',data_type);
    hvol = refreshScreen(hvol);
    
    % Store the data as nifti and check against the segmentation to see if
    % parameter nifti aligns with the segmentation:
    functionals2nifti(hvol,1 , fname);
    mprfCheckParameterNiftiAlignment(cls, fname);
    
    prf_par_exp.(cur_param) = rmGet(hvol.rm.retinotopyModels{1},cur_param);
    
    switch lower(cur_param)
        
        case 'beta'
            
            % Compute the maximum response for every included pRF, by
            % reconstructing the pRF, multiplying it with the stimulus and
            % it's beta, and taking the maximum response from the
            % predicted time series
            mresp = mprfComputeMaximumResponse(rm_stim,sigma_us,x0,y0,prf_par_exp.(cur_param),sm_mask);
                 
            % Store the maximum responses as a nifti file:
            
            fname = fullfile(mrVNif_prf_data,'mresp.nii.gz');
            prf_par_exp.('mresp') = mresp;
            hvol = viewSet(hvol,'map',{mresp});
            functionals2nifti(hvol, 1, fname);
            mprfCheckParameterNiftiAlignment(cls, fname);
            
            prf_par_exp.('mresp') = mresp;
            
             % Smooth the maximum responses on the cortical surface
            [mresp_sm, wConMat] = dhkGraySmooth(hvol,mresp,[ ],wConMat, sm_mask);
            
            % Export smoothed maximum responses as a nifti:
            hvol = viewSet(hvol,'map',{mresp_sm});
            fname = fullfile(mrVNif_prf_data,'mresp_smoothed.nii.gz');
            functionals2nifti(hvol, 1, fname);
            mprfCheckParameterNiftiAlignment(cls, fname);
            
            prf_par_exp.('mresp_smoothed') = mresp_sm;
            
            
            % Recompute the beta by dividing the smoothed maxumimum
            % responses by the maximum response given the stimulus and
            % smoothed pRF parameters:
            recomp_beta = mprfRecomputeBetas(rm_stim,sigma_smooth,x0_smooth,y0_smooth,mresp_sm);
            
            % Store as nifti:
            fname = fullfile(mrVNif_prf_data,'recomp_beta.nii.gz');
            hvol = viewSet(hvol,'map',{recomp_beta});
            functionals2nifti(hvol, 1, fname);
            mprfCheckParameterNiftiAlignment(cls, fname);
            
            prf_par_exp.('recomp_beta') = recomp_beta;
            
            
        case {'x','y','sigma'}
            
            % Smooth the current paramter:
            [tmp_sm_par, wConMat] = dhkGraySmooth(hvol,prf_par_exp.(cur_param),[ ],wConMat, sm_mask);
            
            % Export smoothed data as nifti:
            fname = fullfile(mrVNif_prf_data,[cur_param '_smoothed.nii.gz']);
            
            hvol = viewSet(hvol,'map',{tmp_sm_par});
            functionals2nifti(hvol,1 , fname);
            mprfCheckParameterNiftiAlignment(cls, fname);
            
            prf_par_exp.([cur_param '_smoothed']) = tmp_sm_par;
            
            % Store the current parameters as we need them for the Beta
            % computations above:
            if strcmpi(cur_param,'x')
                x0 = prf_par_exp.(cur_param);
                x0_smooth = tmp_sm_par;
            elseif strcmpi(cur_param,'y')                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
                y0  = prf_par_exp.(cur_param);
                y0_smooth = tmp_sm_par;
            elseif strcmpi(cur_param,'sigma')
                sigma_us = prf_par_exp.(cur_param);
                sigma_smooth = tmp_sm_par;
            end
            
    end
end

fname = fullfile(prf_mat,'exported_prf_params.mat');
save(fname, 'prf_par_exp');

% 
% %% building mrMesh
% 
% % Build mesh
% hvol = meshBuild(hvol,'left'); MSH = meshVisualize(viewGet(hvol,'Mesh')); hvol = viewSet(hvol, 'Mesh', MSH); clear MSH;
% 
% % Smooth the mesh
% hvol = viewSet(hvol, 'Mesh', meshSmooth( viewGet(hvol, 'Mesh'), 1));
% 
% % Load Roi local. for shared change the second 1 to 0
% hvol = loadROI(hvol, 'dialog', [],[],1,1); hvol = selectCurROISlice(hvol); hvol = refreshScreen(hvol,0);
% 
% % Update mesh
% hvol = meshColorOverlay(hvol);



%% Export parameters and ROI from mrVista to freesurfer space

%%%%%%%%%%%%%%%%%%%%%%%%%.
% Exporting pRF parameters 
%%%%%%%%%%%%%%%%%%%%%%%%%.

FS_surface = fprintf('/mnt/storage_2/MEG/Retinotopy/Data/Freesurfer_directory/%s/surf',sub);
fs_prf_data = fprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/prf_data/surface/freesurfer',sub);

mmPerVox = viewGet(hvol,'mmpervox');

% White surfaces were used in original code. But since we are using pial
% surface later when exporting parameters from freesurfer to brainstorm
% surfaces, it might be wise to chose pial surface here. Don't think it
% will make much difference because the vertices are same for both white
% and pial surface. Only the coordinate values will change for ex in the
% rois.
surfaces_to_load = {'lh.pial','rh.pial'};

data_file = dir(fullfile(prf_mat,'*.mat'));
load(fullfile(prf_mat,data_file.name));

% Loop over all the parameters stored in the exported data file:
if ~exist('prf_par_exp', 'var')
    par_names = fieldnames(prfParamsExp);
    prf_par_exp = prfParamsExp;
else
    par_names = fieldnames(prf_par_exp);
end

% Keep track of the X and Y variables, both smoothed and unsmoothed to
% compute the eccenricity and polar angle maps as well:
has_x = false;
has_y = false;
has_x_sm = false;
has_y_sm = false;

for n = 1:length(surfaces_to_load)
    
    cur_surf = fullfile(FS_surface,surfaces_to_load{n});
    
    tmp = strsplit(surfaces_to_load{n},'.');
    cur_hs = tmp{1};
    fprintf('Exporting parameters for %s hemisphere:\n',cur_hs);
    
    % Load mesh using fs_meshFromSurface, this creates a mrVista compatible
    % mesh. Using 'my own' function that skips the smoothing:
    mrv_msh = mprf_fs_meshFromSurface(cur_surf);
    fnum = numel(mrv_msh.triangles);
    
    % compute mapping using mrmMapVerticesToGray (mrmMapGrayToVertices):
    cur_v2gmap = mrmMapVerticesToGray(mrv_msh.vertices, viewGet(hvol,'nodes'),...
        mmPerVox);
    
    % Add the vertexGrayMap field to mesh properties
    mrv_msh = meshSet(mrv_msh,'vertexGrayMap',cur_v2gmap);
    
    good_mapping = cur_v2gmap > 0;
    
    for nn = 1:length(par_names)
        
        cur_par_name = par_names{nn};
        tmp = nan(size(cur_v2gmap));
        
        
        if sum(double(cur_par_name)) == sum(double('beta'))
            tmp_data = squeeze(prf_par_exp.(cur_par_name)(:,:,1));
            tmp(good_mapping) = tmp_data(cur_v2gmap(good_mapping));
            
        else
            tmp(good_mapping) = prf_par_exp.(cur_par_name)(cur_v2gmap(good_mapping));
            
        end
        
        % Check if we have x, y, x_smoothed or y_smoothed:
        if regexp(cur_par_name,'x\>');
            has_x = true;
            x_pos = tmp;
            disp(cur_par_name)
        elseif regexp(cur_par_name,'y\>');
            has_y = true;
            y_pos = tmp;
            disp(cur_par_name)
        elseif regexp(cur_par_name,'x_smoothed\>');
            has_x_sm = true;
            x_pos_sm = tmp;
            disp(cur_par_name)
        elseif regexp(cur_par_name,'y_smoothed\>');
            has_y_sm = true;
            y_pos_sm = tmp;
            disp(cur_par_name)
            
            
        end
        % If we have the necessary data to compute polar angle and
        % eccentricity, do that as well and store the results:
        if has_x && has_y
            [tmp_ang, tmp_ecc] = cart2pol(x_pos, y_pos);
            
            fname = fullfile(fs_prf_data,[cur_hs '.polar_angle']);
            write_curv( fname,tmp_ang, fnum);
            fprintf('Polar angle\n')
            
            fname = fullfile(fs_prf_data,[cur_hs '.eccentricity']);
            write_curv( fname,tmp_ecc, fnum);
            fprintf('Eccentricity\n')
            has_x = false;
            has_y = false;
            
        end
        
        % Same for smoothed versions:
        if has_x_sm && has_y_sm
            [tmp_ang, tmp_ecc] = cart2pol(x_pos_sm, y_pos_sm);
            
            
            fname = fullfile(fs_prf_data,[cur_hs '.polar_angle_smoothed']);
            write_curv( fname,tmp_ang, fnum);
            fprintf('Polar angle smoothed\n')
            
            
            fname = fullfile(fs_prf_data,[cur_hs '.eccentricity_smoothed']);
            write_curv( fname,tmp_ecc, fnum);
            fprintf('Eccentricity smoothed\n')
            
            has_x_sm = false;
            has_y_sm = false;
            
        end
        
        
        % output the results:
        fname = fullfile(fs_prf_data,[cur_hs '.' cur_par_name]);
        write_curv(fname,tmp, fnum);
        fprintf('%s\n', cur_par_name)
    end
    
end

   
%%%%%%%%%%%%%%%%%%%%%%%%%.
% Exporting ROI 
%%%%%%%%%%%%%%%%%%%%%%%%%.
   
% exports ROIs drawn on a Mesh surface to the Freesurfer
% space. Dilates the ROIs a little bit (2 iterations), as they look very 'patchy'
% after the initial mapping to the surface.

% Both hemisphere:
hs = {'LEFT','RIGHT'};
roi_dir = fprintf('/mnt/storage_2/MEG/Retinotopy/Data/fMRI/%s/Gray/ROIs',sub);
fs_roi_dir = fprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/rois/surface/freesurfer',sub);
bs_roi_dir = fprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/rois/surface/brainstorm',sub);

for nn = 1:length(surfaces_to_load)
    if strfind(surfaces_to_load{nn},'lh.')
        roi_path = fullfile(roi_dir);
        tmp = dir(fullfile(roi_path,'L*.mat'));
        roi_names = {tmp.name};
        hs_tag = 100;
    elseif strfind(surfaces_to_load{nn}, 'rh.')
        roi_path = fullfile(roi_dir);
        tmp = dir(fullfile(roi_path,'R*.mat'));
        roi_names = {tmp.name};
        hs_tag = 200;
    end
    
    
    if iscell(roi_names)
        
    else
        roi_names = {roi_names};
        
    end
    
    if iscell(roi_path)
        
    else
        roi_path = {roi_path};
    end
    
    
    % Create a mesh from freesurfer:
    [~,cur_hs] = fileparts(surfaces_to_load{nn});
    cur_surf = fullfile(FS_surface,surfaces_to_load{nn});
    mrv_msh = mprf_fs_meshFromSurface(cur_surf);
    
    % Add the vertexGrayMap field to mesh properties
    %mrv_msh = meshSet(mrv_msh,'vertexGrayMap',cur_v2gmap);
    
    % Add mesh to the volume view:
    hvol = viewSet(hvol,'add and select mesh',mrv_msh);
    
    % Set the roiDilateIterations parameter te prevent 'patchy' outcomes
    prefs = mrmPreferences;
    prefs.roiDilateIterations = 2;
    mrmPreferences(prefs);
    
    % Keep track of the ROI vertices. We de this two times: one for the
    % original ROI vertices, and one for the dilated ROI vertices.
    % Dilation will cause ROIs to claim vertices of adjacent ROIs. We use
    % these variables to keep track of which vertices originally belonged
    % to which ROI and prune the ROIs later on:
    roi_vert_inds_all = cell(1,length(roi_names));
    roi_vert_inds_all2 = cell(1,length(roi_names));
    
    cnt = 0;
    fprintf('Mapping ROIs to freesurfer surface\n')
    % Loop over the ROIS
    for n = 1:length(roi_names)
        [~,r_name] = fileparts(roi_names{n});
        cnt = cnt+1;
        
        % Load the current ROI
        load(fullfile(roi_path{1},roi_names{n}))
        % Get the ROI indices in the gray graph and load the nodes:
        [~, idx] = intersectCols(hvol.coords, ROI.coords);
        roi_nodes = hvol.nodes(:,idx);
        
        % Map the gray nodes to the surface vertices, the results are
        % much cleaner than when mapping the other way around.
        cur_mapping = mrmMapGrayToVertices(roi_nodes, mrv_msh.vertices, mmPerVox, 2);
        
        % Not all voxels can be mapped, give a warning of more than 25% but
        % less than 95% could not be mapped.
        missing = mean(cur_mapping == 0);
        if missing  > .25 && missing < .95
            warning('Could not map %0.2f percent of voxels to surface for %s',...
                missing *100,...
                r_name)
        end
        
        if  missing >= .95
            warning('Could not map %0.2f percent of voxels to surface for %s.\nSKIPPING THIS ROI',...
                missing *100,...
                r_name)
            continue
        end
        
        % ROIs can look patchy on a cortical surface, so dilate them a
        % little be to remove holes in the ROI
        % Keep the orignal and dilated ROI indices
        roiVertInds = unique(cur_mapping(cur_mapping > 0));
        roiVertInds2 = adjustPerimeter(roiVertInds, 0, hvol, prefs);
        
        roi_vert_inds_all{n} = roiVertInds;
        roi_vert_inds_all2{n} = roiVertInds2;
        
    end
    
    % Prune the ROIs to remove overlap due to dilate, keeping the original
    % boundaries between ROIS:
    pruned_roi_inds = mprfResolveDilatedROIOverlap(roi_vert_inds_all2, roi_vert_inds_all);
    fprintf('Done.\n')
    
    all_rois = nan(1,length(mrv_msh.vertices));
    
    fprintf('Combining ROIs\n')
    
    
    cnt = 0 + hs_tag;
    % Now, store the pruned ROI indices, both as a separate file
    for n = 1:length(roi_names)
        cnt = cnt+1;
        roi_tag = mprf__get_roi_tags(roi_names{n});
        
        tag_to_idx.(roi_tag) = cnt;
        
        [~,roi_name] = fileparts(roi_names{n});
        tmp = nan(1,length(mrv_msh.vertices));
        
        if numel(pruned_roi_inds{n})
            
            all_rois(pruned_roi_inds{n}) = cnt;
            tmp(pruned_roi_inds{n}) = cnt;
            
            out_file_name = [cur_hs '.' roi_name(2:end)];
            fname = fullfile(fs_roi_dir,out_file_name);
            
            write_curv(fname, tmp,1);
            
        else
            [~,r_name] = fileparts(roi_names{n});
            warning('No indices found for %s, SKIPPING',...
                r_name);
        end
        
        
    end
    
    fprintf('Done.\n')
    out_file_name = [cur_hs '.all_rois'];
    
    if mean(isnan(all_rois)) == 1
        warning('The combined ROI could not be created for %s',cur_hs)
        
    else
        fname = fullfile(fs_roi_dir,out_file_name);
        write_curv(fname, all_rois,1);
    end
    
    mask = isnan(all_rois);
    out_file_name2 = [cur_hs '.all_rois_mask'];
    
    fname = fullfile(fs_roi_dir,out_file_name2);
    write_curv(fname, mask,1);
    
end

fs_tag_dir = strcat(fs_roi_dir);
save(fullfile(fs_tag_dir, 'fs_tag_to_idx'));

bs_tag_dir = strcat(bs_roi_dir);
save(fullfile(bs_tag_dir, 'bs_tag_to_idx'));


%% Export parameters and ROI from freesurfer to brainstorm space

bs_model_path = fprintf('/mnt/storage_2/MEG/Retinotopy/Data/Brainstorm_db/data/%s/R0774_RETMEG_Block1_5.12.17/',sub); 
bs_anat_path = fprintf('/mnt/storage_2/MEG/Retinotopy/Data/Brainstorm_db/anat/%s/',sub);

bs_model_file = strcat(bs_model_path,'headmodel_surf_os_meg.mat');
bs_anat_file = strcat(bs_anat_path,'subjectimage_T1.mat');

bs_prf_data = fprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/prf_data/surface/brainstorm',sub);
bs_roi_dir = fprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/rois/surface/brainstorm',sub);

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
    
elseif any(strcmpi('white',strsplit(bs_head_model_surf,'_')))
    surfaces_to_load = {'white'};
    
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
        surf_file = fullfile(FS_surface,cur_surf);
        
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
    
    pname = fs_prf_data;
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
                fname = fullfile(bs_prf_data,cur_out_file);
                
            elseif strcmpi(w_export,'roi')
                fname = fullfile(bs_roi_dir,cur_out_file);
                
            end
            
            write_curv(fname,both_bs_data_out,1);
            fprintf('%s\n',cur_out_file);
        end
        
    end    
    
    
    
    %%%%%%%%%%%%%%%%%%
    % Export rois    %
    %%%%%%%%%%%%%%%%%%
    
    pname = fs_roi_dir;
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
                fname = fullfile(bs_prf_data,cur_out_file);
                
            elseif strcmpi(w_export,'roi')
                fname = fullfile(bs_roi_dir,cur_out_file);
                
            end
            
            write_curv(fname,both_bs_data_out,1);
            fprintf('%s\n',cur_out_file);
        end
        
    end    

        
end
        
        
%% Make predictions for every vertex in brainstorm surface
% Make predictions for every brainstorm vertex
% pred_resp = mprf__predicted_prf_response(model, stimulus, prf, roi, iter);
%
% input - model (information about the different options that describes the model)
%       - MEG stimulus
%       - prf parameters on the brainstorm space
%       - roi vertices on the brainstorm surfaces if available
%
% output - predicted response time series for every vertex on brainstorm
%          surface

% Folder to save output (predicted response at brainstorm vertices to MEG
% stimulus
bs_pred_resp = sprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/pred_resp_bs',sub);

% MEG stimulus
% CAN ALSO BE OBTAINED FROM MEG DATA FOLDER
stim_file = sprintf('/mnt/storage_2/MEG/Retinotopy/Subject_sessions/%s/stimuli/meg/imported_stimulus/meg_stimulus.mat',sub);
tmp_stimulus = load(stim_file);
stimulus  = tmp_stimulus.meg_stim;

% brainstorm surface files folder 
bs_surf_data = sprintf('/mnt/storage_2/MEG/Retinotopy/Subject_sessions/%s/prf_data/surface/brainstorm',sub);

% describe the model
model_type = 'original (phase ref amplitude) (model fit) (leave one out)';
model = mprf_generate_model_params(model_type);

% path for the prf parameters on the brainstorm surface
fpath_x0 = fullfile(bs_surf_data,'pial.x_smoothed');
fpath_y0 =  fullfile(bs_surf_data,'pial.y_smoothed');
fpath_sigma = fullfile(bs_surf_data,'pial.sigma_smoothed');
fpath_beta = fullfile(bs_surf_data,'pial.recomp_beta');
fpath_ve = fullfile(bs_surf_data,'pial.varexplained');

% pRF parameters on brainstorm surface
prf.x0.val = read_curv(fpath_x0);
prf.y0.val = read_curv(fpath_y0);
prf.sigma.val = read_curv(fpath_sigma);
prf.beta.val = read_curv(fpath_beta);
prf.ve.val = read_curv(fpath_ve);

% predicted pRF response time series on brainstorm surface

% Not using any roi mask. Should be changes to include ROIs (Wang atlas)
roi.mask = ones(size(prf.x0.val));
roi.idx_out = ones(size(prf.x0.val));

rois = 1;
pred_resp = {zeros(size(stimulus.im,2), size(roi.mask,1),numel(rois))};

nn = 1;% number of rois 
nnn = 1;% number of bootstrapping iterations
cur_in = logical(roi.mask);


if model.params.beta_thr
    if model.params.beta_thr_vals(1) == 0
        range = [prctile(prf.beta.val,0) - 1 prctile(prf.beta.val,model.params.beta_thr_vals(2))];
        
    else
        range = [prctile(prf.beta.val,model.params.beta_thr_vals(1)) prctile(prf.beta.val,model.params.beta_thr_vals(2))];
        
    end
    cur_in = cur_in & prf.beta.val > range(1) & prf.beta.val < range(2);
    
end

if model.params.ve_thr
    range = [model.params.ve_thr_vals(1) model.params.ve_thr_vals(2)];
    cur_in = cur_in & prf.ve.val > range(1) & prf.ve.val < range(2);
end

orig_idx = zeros(size(cur_in));
orig_idx(cur_in) = find(cur_in);
for n = 1:1000:size(prf.sigma.val,1)
    
    end_idx = n+999;
    if end_idx > size(prf.sigma.val,1)
        end_idx = size(prf.sigma.val,1);
    end
    
    
    cur_good_data = cur_in(n:end_idx);
    cur_orig_idx = orig_idx(n:end_idx);
    
    cur_sigma = prf.sigma.val(n:end_idx);
    cur_sigma = cur_sigma(cur_good_data);
    
    
    cur_X0 = prf.x0.val(n:end_idx);
    cur_X0 = cur_X0(cur_good_data);
    
    
    cur_Y0 = prf.y0.val(n:end_idx);
    cur_Y0 = cur_Y0(cur_good_data);
    
    cur_beta = prf.beta.val(n:end_idx);
    cur_beta = cur_beta(cur_good_data);
    
    cur_orig_idx = cur_orig_idx(cur_good_data);
    
    RF = rfGaussian2d(stimulus.X,stimulus.Y,cur_sigma,cur_sigma,false, cur_X0,cur_Y0);
    cur_pred = bsxfun(@times, stimulus.im' * RF,cur_beta');
    
    pred_resp{nnn}(:,cur_orig_idx,nn) = cur_pred;
    
end

disp(sum(sum(isnan(pred_resp{1}))));

% Plot predicted response BS surface
figure, plot(nanmean(pred_resp{1},2))

% save predicted response to MEG stimuli for every brainstorm vertex
save(fullfile(bs_pred_resp,'/pred_resp_bs'),'pred_resp');

% pred_resp will now be converted to every meg sensor by multiplying it
% with the gain matrix

%% Predicted response at MEG sensors to MEG stimuli
% Multiply with gain matrix to get the predictions per vertex
% meg_resp = mprf__compute_predicted_meg_response(bs, pred_resp, channels);
% input - gain matrix
%       - predicted response time series for every vertex on brainstorm
%         surface 
%       - channel information 
%
% output - predicted response time series for every channel on MEG surface

% head model
% CAN ALSO BE OBTAINED FROM BRAINSTORM FOLDER DIRECTLY
bs.model_file = sprintf('/mnt/storage_2/MEG/Retinotopy/Subject_sessions/%s/source/brainstorm/head_model/head_model_tess_cortex_pial_low.mat',sub);
bs_model = load(bs.model_file);
bs.lead_field = bst_gain_orient(bs_model.Gain,bs_model.GridOrient);
bs.lead_field2 = bs.lead_field(~isnan(bs.lead_field(:,1)),:); % Remove NaNs....
bs.keep_sensors = ~isnan(bs.lead_field(:,1));

% predicted response at brainstorm vertices to MEG stimulus
bs_pred_resp = sprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/pred_resp_bs',sub);
load(fullfile(bs_pred_resp,'pred_resp_bs.mat'));

% MEG  Channels
channels.data = 0:156; % MEG channels in which MEG data is acquired
channels.trigger = 160:167; % MEG channels used as triggers
channels.diode = 191; % Channel for collecting diode values



if iscell(pred_resp)
    sz_cell = size(pred_resp,2);
    meg_resp = cell(size(pred_resp));
    
    for n = 1:sz_cell
        fprintf('Iteration: %d\n',n)
        cur_pred = pred_resp{n};
        
        
        if size(cur_pred,3) == 1
            meg_resp{n} = cur_pred * bs.lead_field(bs.keep_sensors,:)';
            meg_resp{n} = meg_resp{n}(:,channels.data+1);
            
        else
            % loop over ROIs here
            error('Not implemented')
        end
        
    end
end   



%%

meg_data_file_path = sprintf('/mnt/storage_2/MEG/Retinotopy/Subject_sessions/%s/data/meg/preproc/pp/epoched_data_hp_preproc_denoised.mat',subjid);
plot_figure = 1;

pred.pred_resp = pred_resp;
pred.prf = prf;
pred.bs = bs;
pred.roi = roi;
pred.model = model;
pred.stimulus = stimulus;
pred.syn = 0;
pred.meg_resp = meg_resp;
pred.channels = channels;
cur_time = datestr(now);
cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
pred.cur_time = cur_time;

mprfSession_run_original_model_server(pred,meg_data_file_path,plot_figure);

