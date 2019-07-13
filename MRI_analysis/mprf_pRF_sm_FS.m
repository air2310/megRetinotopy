function mprf_pRF_sm_FS(dirPth)



%% ----------
% File paths
% ----------
anat_dir = fullfile(dirPth.fmri.mrvPth, '3DAnatomy'); 
anat_file = fullfile(anat_dir,'t1.nii.gz');
roi_dir  = strcat(anat_dir,'/ROIs/');

freesurfer_surface = dirPth.fs.surfPth;

mrSession_dir = dirPth.fmri.mrvPth; 

prf_data_mrVmat = strcat(dirPth.fmri.saveDataPth_prfMrv,'/mat');

% directories to save results
% original and smoothed pRF parameters in mrVista space in .nii 
prf_dir_FS = dirPth.fmri.saveDataPth_prfFS;
roi_dir_FS = dirPth.fmri.saveDataPth_roiFS;
roi_dir_BS = dirPth.fmri.saveDataPth_roiBS;
% ----------


% step inside the vistasession directory contain mrSESSION.mat
cd(mrSession_dir);
% We need a volume view:
data_type = 'Averages';
setVAnatomyPath(anat_file);
hvol = initHiddenGray;
% Set the volume view to the current data type and add the RM model
hvol = viewSet(hvol,'curdt',data_type);

% Load mrVista retinotopy Gray file 
rm_model = dirPth.fmri.vistaGrayFitFile;
hvol = rmSelect(hvol,1,rm_model);

mmPerVox = viewGet(hvol,'mmpervox');

% White surfaces were used in original code. But since we are using pial
% surface later when exporting parameters from freesurfer to brainstorm
% surfaces, it might be wise to chose pial surface here. Don't think it
% will make much difference because the vertices are same for both white
% and pial surface. Only the coordinate values will change for ex in the
% rois.
surfaces_to_load = {'lh.pial','rh.pial'};

data_file = dir(fullfile(prf_data_mrVmat,'*.mat'));
load(fullfile(prf_data_mrVmat,data_file.name));

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
    
    cur_surf = fullfile(freesurfer_surface,surfaces_to_load{n});
    
    tmp = strsplit(surfaces_to_load{n},'.');
    cur_hs = tmp{1};
    fprintf('Exporting parameters for %s hemisphere:\n',cur_hs);
    
    % Load mesh using fs_meshFromSurface, this creates a mrVista compatible
    % mesh. Using 'my own' function that skips the smoothing:
    mrv_msh = mprf_fs_meshFromSurface(cur_surf);
    fnum = numel(mrv_msh.triangles);
    
    % compute mapping using mrmMapVerticesToGray (mrmMapGrayToVertices):
    cur_v2gmap = mrmMapVerticesToGray(mrv_msh.vertices, viewGet(hvol,'nodes'),...
        mmPerVox,[],5);
    
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
        if strcmpi(cur_par_name,'x')
            has_x = true;
            x_pos = tmp;
            disp(cur_par_name)
        elseif strcmpi(cur_par_name,'y')
            has_y = true;
            y_pos = tmp;
            disp(cur_par_name)
        elseif strcmpi(cur_par_name,'x_smoothed')
            has_x_sm = true;
            x_pos_sm = tmp;
            disp(cur_par_name)
        elseif strcmpi(cur_par_name,'y_smoothed')
            has_y_sm = true;
            y_pos_sm = tmp;
            disp(cur_par_name)
            
            
        end
        % If we have the necessary data to compute polar angle and
        % eccentricity, do that as well and store the results:
        if has_x && has_y
            [tmp_ang, tmp_ecc] = cart2pol(x_pos, y_pos);
            
            fname = fullfile(prf_dir_FS,[cur_hs '.polar_angle']);
            write_curv( fname,tmp_ang, fnum);
            fprintf('Polar angle\n')
            
            fname = fullfile(prf_dir_FS,[cur_hs '.eccentricity']);
            write_curv( fname,tmp_ecc, fnum);
            fprintf('Eccentricity\n')
            has_x = false;
            has_y = false;
            
        end
        
        % Same for smoothed versions:
        if has_x_sm && has_y_sm
            [tmp_ang, tmp_ecc] = cart2pol(x_pos_sm, y_pos_sm);
            
            
            fname = fullfile(prf_dir_FS,[cur_hs '.polar_angle_smoothed']);
            write_curv( fname,tmp_ang, fnum);
            fprintf('Polar angle smoothed\n')
            
            
            fname = fullfile(prf_dir_FS,[cur_hs '.eccentricity_smoothed']);
            write_curv( fname,tmp_ecc, fnum);
            fprintf('Eccentricity smoothed\n')
            
            has_x_sm = false;
            has_y_sm = false;
            
        end
        
        
        % output the results:
        fname = fullfile(prf_dir_FS,[cur_hs '.' cur_par_name]);
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
    cur_surf = fullfile(freesurfer_surface,surfaces_to_load{nn});
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
            fname = fullfile(roi_dir_FS,out_file_name);
            
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
        fname = fullfile(roi_dir_FS,out_file_name);
        write_curv(fname, all_rois,1);
    end
    
    mask = isnan(all_rois);
    out_file_name2 = [cur_hs '.all_rois_mask'];
    
    fname = fullfile(roi_dir_FS,out_file_name2);
    write_curv(fname, mask,1);
    
end

fs_tag_dir = strcat(roi_dir_FS);
save(fullfile(fs_tag_dir, 'fs_tag_to_idx'));

bs_tag_dir = strcat(roi_dir_BS);
save(fullfile(bs_tag_dir, 'bs_tag_to_idx'));