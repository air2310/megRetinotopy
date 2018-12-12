function megRet_FS2BS(s, roiType)

if ~exist('MRIread')
    addpath(genpath('/Applications/freesurfer/matlab'));
    addpath(genpath('/Applications/freesurfer/fsfast/toolbox'));
end

% open brainstorm
brainstorm

% Load the brainstorm head model file, as we need to know on which surface
% it is defined:
% bsOverlay = load(fullfile(s.BS.surface.pth, 'tess_cortex_pial_low.mat'));

% Set the freesurfer subject dir
fsSubjDir = fullfile(s.freeSurferDir, s.fsSubject);
setenv('SUBJECTS_DIR', s.freeSurferDir)

% Set the brainstorm subject dir
bsSubjDir = s.BS.surface.pth;

% Save ROIs
saveRoiDir = fullfile(s.outPut.pth, 'rois');
if ~exist(saveRoiDir,'dir'); mkdir(saveRoiDir); end;

% Get the Wang Atlas
prfParamfname = @(hem, type)(sprintf('%s/surf/%s.wang2015_atlas.mgz', fsSubjDir, hem));
surfdat         = @(lh, rh, rc)(setfield(setfield([], 'lh', lh.vol(:)), 'rh', rc*rh.vol(:)));

subFSWang       = surfdat(MRIread(prfParamfname('lh')), MRIread(prfParamfname('rh')), 1);

switch roiType
    
    case 'allRoisWangAtlas'
        
        % get all rois in one data structure
        lh.ROI_data = subFSWang.lh;
        rh.ROI_data = subFSWang.rh;
        
        % get all rois as a mask (0's and 1's)
        lh.ROI_data_mask = subFSWang.lh>0;
        rh.ROI_data_mask = subFSWang.rh>0;

        
        % do interpolation
        wangROI_BS = tess_fs2bst(s.bsSubject, fsSubjDir, lh.ROI_data, rh.ROI_data);
        wangROI_BS_mask = double(tess_fs2bst(s.bsSubject, fsSubjDir, lh.ROI_data_mask, rh.ROI_data_mask));
        wangROI_BS_mask(wangROI_BS_mask==0)=NaN;
        
        % Create file names for interpolated data
        subBSAllROImgzfile = sprintf('%s/pial.all_rois.mgz', saveRoiDir);
        subBSAllROIMaskmgzfile = sprintf('%s/pial.all_rois_mask.mgz', saveRoiDir);

        subBSAllROImatfile = sprintf('%s/allWangROIs_overlay.mat', saveRoiDir);
        
        % save ROIs as mgz file
        MRIwrite(struct('vol', wangROI_BS), subBSAllROImgzfile);
        MRIwrite(struct('vol', wangROI_BS_mask), subBSAllROIMaskmgzfile);
        
        % save ROIs as BS binary files
        cur_dir = pwd; 
        if ~exist(fullfile(saveRoiDir, 'surface', 'brainstorm'), 'dir'); mkdir(fullfile(saveRoiDir, 'surface', 'brainstorm')); end;
        
        cd(fullfile(saveRoiDir, 'surface', 'brainstorm'))

        write_curv('pial.all_rois',  wangROI_BS, 1)
        write_curv('pial.all_rois_mask', wangROI_BS_mask, 1)
        cd(cur_dir)
        
        % save ROIs as mat file
        save(subBSAllROImatfile, 'wangROI_BS');
         
        
        
    case 'individual'
        % in case of individual rois???
        
end

%% Separate PRF data into Left and Right Hemi, then project to FS mesh?

% location of smoothed prf data
loadPrfDataDir = fullfile(s.outPut.pth, 'prf_data');
savePrfDataDirFS = fullfile(loadPrfDataDir, 'surface', 'freesurfer');
savePrfDataDirBS = fullfile(loadPrfDataDir, 'surface', 'brainstorm');


load(fullfile(loadPrfDataDir, 'data_dir', 'exported_prf_params.mat'), 'prf_par_exp');

% Get all the parameters stored in the prf parameter data file:
parNames = fieldnames(prf_par_exp);

% Keep track of the X and Y variables, both smoothed and unsmoothed to
% compute the eccenricity and polar angle maps as well:
has_x = false;
has_y = false;
has_x_sm = false;
has_y_sm = false;

% surfaces to transform from Vista to FS
fsSurfaces = {'lh.pial','rh.pial'};
% fsSurfaces = {'lh.white','rh.white'};


currPath = pwd;
cd(s.vistaSession.pth)
hvol = initHiddenGray;
mmPerVox = viewGet(hvol,'mmpervox');

% Loop over the FS surfaces:
for n = 1:length(fsSurfaces)
    
    currentSurface = fullfile(s.freeSurferDir, s.fsSubject, 'surf', fsSurfaces{n});
    tmp = strsplit(fsSurfaces{n},'.');
    currentHemi = tmp{1};
    fprintf('Exporting parameters for %s hemisphere:\n',currentHemi)
    
    % Load mesh using fs_meshFromSurface, this creates a mrVista compatible
    % mesh. Using 'my own' function that skips the smoothing:
    mrVistaMesh = mprf_fs_meshFromSurface(currentSurface);
    fnum = numel(mrVistaMesh.triangles);
    
    % compute mapping using mrmMapVerticesToGray (mrmMapGrayToVertices):
    currentVista2Graymap = mrmMapVerticesToGray(mrVistaMesh.vertices, viewGet(hvol,'nodes'),...
        mmPerVox);
    
    goodVerticesForMapping = currentVista2Graymap > 0;
    
    % cur_v2gmap is really just an index variable into the gray nodes
    % that tells which gray node is closest to which mesh node. As it uses
    % nearpoints (i.e. nearest neighbour interpolation), I suspect this only
    % works for mrVista based meshes.
    % Would be nice to do something more intelligent than nearest neighbour
    % interpolation, like checking other layers of the Gray graph as well.
    
    % And use the cur_v2gmap to select the pRF parameters of the Gray nodes
    % that are closest to the mesh nodes/vertices. However, cur_v2gmap contains
    % zeros, so we need to deal with them:
    
    % Loop over the parameters:
    for nn = 1:length(parNames)
        
        currentParameterName = parNames{nn};
        tmp = nan(size(currentVista2Graymap));
        
        if sum(double(currentParameterName)) == sum(double('beta'))
            tmp_data = squeeze(prf_par_exp.(currentParameterName)(:,:,1));
            tmp(goodVerticesForMapping) = tmp_data(currentVista2Graymap(goodVerticesForMapping));
            
        else
            tmp(goodVerticesForMapping) = prf_par_exp.(currentParameterName)(currentVista2Graymap(goodVerticesForMapping));
            
        end
        % Check if we have x, y, x_smoothed or y_smoothed:
        if regexp(currentParameterName,'x\>');
            has_x = true;
            x_pos = tmp;
            disp(currentParameterName)
        elseif regexp(currentParameterName,'y\>');
            has_y = true;
            y_pos = tmp;
            disp(currentParameterName)
        elseif regexp(currentParameterName,'x_smoothed\>');
            has_x_sm = true;
            x_pos_sm = tmp;
            disp(currentParameterName)
        elseif regexp(currentParameterName,'y_smoothed\>');
            has_y_sm = true;
            y_pos_sm = tmp;
            disp(currentParameterName)
            
            
        end
        % If we have the necessary data to compute polar angle and
        % eccentricity, to that as well and store the results:
        if has_x && has_y
            [tmp_ang, tmp_ecc] = cart2pol(x_pos, y_pos);
            
            fname = fullfile(savePrfDataDirFS,[currentHemi '.polar_angle']);
            write_curv( fname,tmp_ang, fnum);
            fprintf('Polar angle\n')
            
            fname = fullfile(savePrfDataDirFS,[currentHemi '.eccentricity']);
            write_curv( fname,tmp_ecc, fnum);
            fprintf('Eccentricity\n')
            has_x = false;
            has_y = false;
            
        end
        
        % Same for smoothed versions:
        if has_x_sm && has_y_sm
            [tmp_ang, tmp_ecc] = cart2pol(x_pos_sm, y_pos_sm);
            
            
            fname = fullfile(savePrfDataDirFS,[currentHemi '.polar_angle_smoothed']);
            write_curv( fname,tmp_ang, fnum);
            fprintf('Polar angle smoothed\n')
            
            
            fname = fullfile(savePrfDataDirFS,[currentHemi '.eccentricity_smoothed']);
            write_curv( fname,tmp_ecc, fnum);
            fprintf('Eccentricity smoothed\n')
            
            has_x_sm = false;
            has_y_sm = false;
            
        end
        % output the results:
        fname = fullfile(savePrfDataDirFS,[currentHemi '.' currentParameterName]);
        write_curv(fname,tmp, fnum);
        fprintf('%s\n', currentParameterName)
        
    end
    
end

cd(currPath)

%% Project to FS data back to brainstorm..
% (we could've gone straight from the smoothing into a mgz file. not sure why this is done this way...)

prfParamsToUse = {'beta', 'eccentricity', 'eccentricity_smoothed', 'mask', ...
    'mresp', 'mresp_smoothed', 'polar_angle', 'polar_angle_smoothed', 'recomp_beta', ...
    'sigma', 'sigma_smoothed', 'varexplained', 'x', 'x_smoothed', 'y', 'y_smoothed'};

% Loop over the param files listed
for nn = 1:length(prfParamsToUse)
    
    thisParam = prfParamsToUse{nn};
    % Get the prf surface
    prfParamfname = @(hem)(sprintf('%s/%s.%s', savePrfDataDirFS, hem, thisParam));
    
    surfdataFS.lh = read_curv(prfParamfname('lh')); 
    surfdataFS.rh = read_curv(prfParamfname('rh'));


    % do interpolation
    surfdataBS = tess_fs2bst(s.bsSubject, fsSubjDir, surfdataFS.lh, surfdataFS.rh);
        
    % Create file names for interpolated data
    surfdataBSfile = sprintf('%s/%s.%s', savePrfDataDirBS, 'pial', thisParam);
    if ~exist(savePrfDataDirFS, 'dir'); mkdir(savePrfDataDirFS); end;
    
    % save ROIs as mgz file
    write_curv(surfdataBSfile,surfdataBS,1);
end  

    
    




















