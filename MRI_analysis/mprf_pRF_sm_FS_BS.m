function mprf_pRF_sm_FS_BS(subjID, dirPth,opt)
% mprf_pRF_sm_FS_BS(subjID, dirPth,plot_stim)
%
% Function to export prf parameters (unsmooth/smooth) and ROIs from
% FreeSurfer (FS) left and right hemisphere surfaces to Brainstorm pial
% surface. Currently, the Brainstorm surface is a downsampled FreeSurfer
% pial surface, with 15002 vertices. It will also save a file where you
% combine the left and hemisphere of the  FreeSurfer surface (note that
% these data are on the midgray surface, not exported to the pial surface
% -- ek: not sure if that matters?)
%
% INPUTS:
%   subjID      :   name of subject, e.g. 'wlsubj004'. (string)
%   dirPth      :   paths locating subject's data and files (struct, see loadPaths.m)
%   opt         :   struct with boolean flags, needed to request using rois
%                   defined in mrVista space

%% ----------
% File paths
% -----------
prf_dir_FS = dirPth.fmri.saveDataPth_prfFS;
roi_dir_FS = dirPth.fmri.saveDataPth_roiFS;
prf_dir_BS = dirPth.fmri.saveDataPth_prfBS;
roi_dir_BS = dirPth.fmri.saveDataPth_roiBS;

if (~exist(prf_dir_BS, 'dir') || ~exist(roi_dir_BS, 'dir') || ~exist(roi_dir_FS, 'dir'))
    mkdir(prf_dir_BS);
    mkdir(roi_dir_BS);
    mkdir(roi_dir_FS);
end
% ----------

%% ------------------------------------------------------------------------
% Exporting pRF parameters from freesurfer vertices to brainstorm vertices
% -------------------------------------------------------------------------

% Open up brainstorm if database can't be found
if ~exist('GlobalData.DataBase.iProtocol', 'var')
    brainstorm nogui;
end

% Parameters and ROIs will be mapped to a downsampled version of the
% Freesurfer pial surface
surfaceToMapTo    = 'pial';

% Get left file names, do we can extract the prf parameter name (could also
% be done with the right file names..)
lhFiles = dir(fullfile(prf_dir_FS,'lh.*'));

for nn = 1:length(lhFiles)
    
    % Get prf parameter name
    curLHFile = lhFiles(nn).name;
    [~,~,paramName] = fileparts(curLHFile);
    paramName = strsplit(paramName, '.');
    paramName = paramName{2};
    
    prfParamfname = @(hem)(sprintf('%s/%s.%s', prf_dir_FS, hem, paramName));
    
    % Load Freesurfer data
    surfdataFS.lh = read_curv(prfParamfname('lh'));
    surfdataFS.rh = read_curv(prfParamfname('rh'));
        
    % load and concatenate:
    bothHemiFSData = [surfdataFS.lh; surfdataFS.rh];
    
    % Find corresponding BS vertices for FS surface data
    bothHemiBSData = tess_fs2bst(subjID, dirPth.fs.segPth, surfdataFS.lh, surfdataFS.rh);
    
    % Save the BS results:
    curFileToSave = [surfaceToMapTo '.' paramName];
    fname = fullfile(prf_dir_BS,curFileToSave);
    write_curv(fname,bothHemiBSData,1);
    fprintf('(%s): Brainstorm combined hemi files: %s\n',mfilename, curFileToSave);
    
    % Save the FS results:
    curFileToSave = [surfaceToMapTo '.' paramName];
    fname = fullfile(prf_dir_FS,curFileToSave);
    write_curv(fname,bothHemiFSData,1);
    fprintf('(%s): Freesurfer combined hemi files: %s\n',mfilename, curFileToSave);
    
end


%% --------------------------------------------------------------------
% Export rois
%----------------------------------------------------------------------
if opt.roimrvToFS % in case you want to use the hand-drawn rois from mrVista
    
    lhFiles = dir(fullfile(roi_dir_FS,'lh.*'));
    
    % Load ROIs drawn in mrVista surface and exported to freesurfer space
    for ii  = 1:length(lhFiles)
        
        paramName = lhFiles{ii};
        paramName = paramName(2:end);
        
        % Find the corresponding rh file:
        curLHFile = ['lh.' paramName];
        curRHFile = ['rh.' paramName];
        
        % load and concatenate:
        surfROIFS.lh = read_curv(fullfile(roi_dir_FS,curLHFile));
        surfROIFS.rh = read_curv(fullfile(roi_dir_FS,curRHFile));
        
        % load and concatenate:
        bothHemiFSROI = [surfROIFS.lh; surfROIFS.rh];
    
        % Find BS vertices
        bothHemiBSROI = tess_fs2bst(subjID, dirPth.fs.segPth, surfROIFS.lh, surfROIFS.rh);

        % Store the results:
        curFileToSave = [surfaceToMapTo '.' paramName];
        fname = fullfile(roi_dir_BS,curFileToSave);
        write_curv(fname,bothHemiBSROI,1);
        fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, curFileToSave);
        
        
        % Create roi mask FS pial file
        fname = fullfile(prf_dir_FS,curFileToSave);
        write_curv(fname,bothHemiFSROI,1);
        fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, curFileToSave);
    end
    
    
else  % in case you want to use Wang et al. (2015) atlas rois
    
    paramName = 'wang2015_atlas.mgz';
    
    roifname = @(hem)(sprintf('%s/%s.%s', dirPth.fs.surfPth, hem, paramName));
    
    tmp = MRIread(roifname('lh'));
    surfROIFS.lh = squeeze(tmp.vol);
    tmp = MRIread(roifname('rh'));
    surfROIFS.rh = squeeze(tmp.vol);
    
    % load and concatenate:
    bothHemiFSROI = [surfROIFS.lh; surfROIFS.rh];
    
    % Find BS vertices
    bothHemiBSROI = tess_fs2bst(subjID, dirPth.fs.segPth, surfROIFS.lh, surfROIFS.rh);
    
    % ---- Store full Wang atlas ----
    % BRAINSTORM Pial file
    curFileToSave = [surfaceToMapTo '.wang2015_atlas'];
    write_curv(fullfile(roi_dir_BS,curFileToSave),bothHemiBSROI,1);
    write_curv(fullfile(prf_dir_BS,curFileToSave),bothHemiBSROI,1);
    fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, curFileToSave);
    
    % COMBINE HEMI FREESURFER file
    write_curv(fullfile(roi_dir_FS,curFileToSave),bothHemiFSROI,1);
    write_curv(fullfile(prf_dir_FS,curFileToSave),bothHemiFSROI,1);
    fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, curFileToSave);
    
    
    % ---- Store mask of Wang atlas ----
    % BRAINSTORM Pial file
    mask = bothHemiBSROI>0;
    curFileToSave = [surfaceToMapTo '.mask'];
    write_curv(fullfile(roi_dir_BS,curFileToSave),mask,1);
    write_curv(fullfile(prf_dir_BS,curFileToSave),mask,1);
    fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, curFileToSave);
    
    % COMBINE HEMI FREESURFER file
    mask = bothHemiFSROI>0;
    curFileToSave = [surfaceToMapTo '.mask'];
    write_curv(fullfile(roi_dir_FS,curFileToSave),mask,1);
    write_curv(fullfile(prf_dir_FS,curFileToSave),mask,1);
    fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, curFileToSave);
    
    % ---- Store just an roi mask of V1-V3 areas in Wang atlas ----
    V123idx = 1:6; % first six rois are "V1v" "V1d" "V2v" "V2d" "V3v" "V3d"
    
    % BRAINSTORM Pial file
    BS_V123mask = ismember(bothHemiBSROI, V123idx);
    curFileToSave = [surfaceToMapTo '.V123mask'];
    write_curv(fullfile(roi_dir_BS,curFileToSave),BS_V123mask,1);
    write_curv(fullfile(prf_dir_BS,curFileToSave),BS_V123mask,1);
    fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, curFileToSave);
    
    % COMBINE HEMI FREESURFER file
    FS_V123mask = ismember(bothHemiFSROI, V123idx);
    write_curv(fullfile(prf_dir_FS,curFileToSave),FS_V123mask,1);
    write_curv(fullfile(roi_dir_FS,curFileToSave),FS_V123mask,1);
    fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, curFileToSave);
    
end

