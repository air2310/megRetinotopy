function mprfSessionExportVistaROIToFreesurferSurface
% This function exports ROIs drawn on a Mesh surface to the Freesurfer
% space. It uses the same mapping between mrVista Gray and the surface as
% the mprfSessionExportPRFDataToFreeSurferSurface.m. Additionally, it
% dilates the ROIs a little bit (2 iterations), as they look very 'patchy'
% after the initial mapping to the surface. I do not know why this is yet.
%
%
%% TO DO:
% Remove file name dependence at line 108 and 140
% Make sure this file also works with ROIs defined on both Hemispheres
%%
% Load mprfSESSION and do the necessary checks:
load(fullfile(pwd,'mprfSESSION.mat'))
global mprfSESSION

cd(mprfSESSION.orig.vista_dir);

if exist('mprf_fs_meshFromSurface','file')
    
else
    addpath(mprfSESSION.root.root_path);
    if exist('mprf_fs_meshFromSurface','file')
        
    else
        error('Could not add root directory to path')
        
    end
end

if exist('mrmMapGrayToVertices','file');
    
else
    try
        tbUse('vistasoft');
    catch ME
        error('Vistasoft is not on your Matlab path')
        
    end
    
end


if ~isfield(mprfSESSION,'orig') || ~isfield(mprfSESSION.orig,'vista_dir')
    error('Cannot find mrVista directory path, mprfSession is not properly initialized')
end



if ~isfield(mprfSESSION,'source') || ~isfield(mprfSESSION.source,'fs_surface_dir') || ...
        ~exist(mprfSESSION.source.fs_surface_dir,'dir') || isempty(dir(mprfSESSION.source.fs_surface_dir));
    
    error('No imported freesurfer surfaces found,  mprfSession is not properly initialized');
end

% Initialize a Gray view:
hvol = initHiddenGray;
mmPerVox = viewGet(hvol,'mmpervox');

% Freesurfer surfaces to use, again only the white surface:
surfaces_to_load = {'lh.white','rh.white'};
% Both hemisphere:
hs = {'LEFT','RIGHT'};


% Loop over the surfaces:
for nn = 1:length(surfaces_to_load)
    
    % Ask for the ROIs we want to export for the current hemisphere:
    [roi_names, roi_path] = uigetfile('*',...
        sprintf('Please select the %s hemisphere ROIs you want to export',hs{nn}),...
        'MultiSelect','on');
    
    if iscell(roi_names)
        
    else
        roi_names = {roi_names};
        
    end
    
    
    if iscell(roi_path)
        
    else
        roi_path = {roi_path};
        
    end
    
    
    
    if ~ischar(roi_names{1});
        fprintf('No valid ROI names provided.\nSkipping %s hemisphere\n',hs{nn})
        
    end
    
    % Create a mesh from freesurfer:
    [~,cur_hs] = fileparts(surfaces_to_load{nn});
    cur_surf = fullfile(mprfSESSION.source.fs_surface_dir, surfaces_to_load{nn});
    mrv_msh = mprf_fs_meshFromSurface(cur_surf);
    
    % Add mesh to the volume view:
    hvol = viewSet(hvol,'add and select mesh',mrv_msh);
    
    %     clrs = mrv_msh.colors;
    %    clrs2 = mrv_msh.colors;
    
    %    surf_data = zeros(1,length(clrs));
    
    %     cm_idx = round(linspace(1,256,length(roi_idx)));
    %     cm = 255 .* colorcube(256);
    
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
        
        fprintf('%s\n',r_name);
        % Keep the orignal and dilated ROI indices
        roiVertInds = unique(cur_mapping(cur_mapping > 0));
        roiVertInds2 = adjustPerimeter(roiVertInds, 0, hvol, prefs);
        
        roi_vert_inds_all{n} = roiVertInds;
        roi_vert_inds_all2{n} = roiVertInds2;
        
    end
    
    fprintf('Done.\n')
    
    fprintf('Pruning ROIs\n')
    
    % Prune the ROIs to remove overlap:
    pruned_roi_inds = mprfResolveDilatedROIOverlap(roi_vert_inds_all2, roi_vert_inds_all);
    fprintf('Done.\n')
    
    all_rois = nan(1,length(mrv_msh.vertices));
    
    fprintf('Combining ROIs\n')
    
    
    cnt = 0;
    % Now, store the pruned ROI indices, both as a separate file
    for n = 1:length(roi_names)
        
        cnt = cnt+1;
        [~,roi_name] = fileparts(roi_names{n});
        tmp = nan(1,length(mrv_msh.vertices));
        
        if numel(pruned_roi_inds{n})
            
            all_rois(pruned_roi_inds{n}) = cnt;
            tmp(pruned_roi_inds{n}) = cnt;
            
            out_file_name = [cur_hs '.' roi_name(2:end)];
            fname = mprfExportDataPath('roi_fs_surface',out_file_name);
            
            write_curv(fname, tmp,1);
            
        else
            [~,r_name] = fileparts(roi_names{n});
            warning('No indices found for %s, SKIPPING',...
                r_name);
        end
        
        
    end
    %         tmp = cm(cm_idx(n),:)';
    %         clrs2(1:3,pruned_roi_inds{n}) = tmp(:,ones(1,length(pruned_roi_inds{n})));
    
    fprintf('Done.\n')
    out_file_name = [cur_hs '.all_rois'];
    
    if mean(isnan(all_rois)) == 1
        warning('The combined ROI could not be created for %s',cur_hs)
        
    else
        fname = mprfExportDataPath('roi_fs_surface',out_file_name);
        write_curv(fname, all_rois,1);
    end
    mask = isnan(all_rois);
    out_file_name2 = [cur_hs '.all_rois_mask'];
    
    fname = mprfExportDataPath('roi_fs_surface',out_file_name2);
    write_curv(fname, mask,1);
    
    
    
    
end




cd(mprfSESSION.init.main_dir);
save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
mrvCleanWorkspace;
clear mprfSESSION

end




