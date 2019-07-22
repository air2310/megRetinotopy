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

hemiToLoad      = {'lh','rh'};        % In this order as BS concatenates the hemispheres like this
surfaceToMap    = 'pial';

fprintf('Exporting surface data:\n');



% Get brainstorm mesh
fs_vertices = [];
% Loop over all (2, left and right) surfaces:
for n = 1:length(hemiToLoad)
    cur_surf = [hemiToLoad{n} '.' surfaceToMap];
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
    cur_out_file = [surfaceToMap '.' par_name];
    fname = fullfile(prf_dir_BS,cur_out_file);
    write_curv(fname,both_bs_data_out,1);
    fprintf('(%s): Brainstorm combined hemi files: %s\n',mfilename, cur_out_file);
    
    % Store the FS results:
    cur_out_file = [surfaceToMap '.' par_name];
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
    
    % Load ROIs drawn in mrVista surface and exported to freesurfer space
    for ii  = 1:length(lh_files)
        par_name = lh_files{ii};
        
        par_name = par_name(2:end);
        % Find the corresponding rh file:
        cur_rh_file = ['rh.' par_name];
        
        % load and concatenate:
        both_fs_data = [read_curv(fullfile(pname,cur_lh_file));...
            read_curv(fullfile(pname,cur_rh_file))];
        
        % Select to correct parameters:
        both_bs_data_out = both_fs_data(bs_vert_idx);
        
        % Store the results:
        cur_out_file = [surfaceToMap '.' par_name];
        fname = fullfile(roi_dir_BS,cur_out_file);
        write_curv(fname,both_bs_data_out,1);
        fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, cur_out_file);
        
        
        % Create roi mask FS pial file
        fname = fullfile(prf_dir_FS,cur_out_file);
        write_curv(fname,both_fs_data,1);
        fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, cur_out_file);
    end
    
    
else
    
    % load wang rois
    combineHemi = true;
    combineRois = false;
    roiData = loadWangROIs(dirPth.fs.surfPth, 'lh.wang2015_atlas*', combineHemi, combineRois);
    
    % Loop over the rois to get separate prf data
    fnROI = fieldnames(roiData);
    for n = 1:length(fieldnames(roiData))
        
        both_fs_data = roiData.(fnROI{n});
        
        % Resample to nearest Brainstorm vertices
        both_bs_data_out = both_fs_data(bs_vert_idx);
        
        % Store full Wang atlas as BS file
        cur_out_file = [surfaceToMap '.' fnROI{n}];
        fname = fullfile(roi_dir_BS,cur_out_file);
        write_curv(fname,both_bs_data_out,1);
        fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, cur_out_file);
        
        % Store combined Wang atlas as FS file
        fname = fullfile(roi_dir_FS,cur_out_file);
        write_curv(fname,both_fs_data,1);
        fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, cur_out_file);
    end
end

