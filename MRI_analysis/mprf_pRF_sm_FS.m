function mprf_pRF_sm_FS(dirPth,opt)
% mprf_pRF_sm_FS(dirPth,plot_stim)
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
% ----------
anatDir = fullfile(dirPth.fmri.mrvPth, '3DAnatomy');
t1file  = fullfile(anatDir,'t1.nii.gz');
roiDir  = strcat(anatDir,'/ROIs/');

fsSurfDir = dirPth.fs.surfPth;

mrSessionDir = dirPth.fmri.mrvPth;

prfDatamrVistaMatDir = strcat(dirPth.fmri.saveDataPth_prfMrv,'/mat');

% directories to save results
% original and smoothed pRF parameters in mrVista space in .nii
prfFSDir = dirPth.fmri.saveDataPth_prfFS;
roiFSDir = dirPth.fmri.saveDataPth_roiFS;
roiBSDir = dirPth.fmri.saveDataPth_roiBS;

if ~exist(prfFSDir, 'dir')
    mkdir(prfFSDir);
    mkdir(roiFSDir);
    mkdir(roiBSDir);
end
% ----------

%% ------------------------------------------------------------------------
% Exporting prf parameters to freesurfer vertices
%--------------------------------------------------------------------------

% step inside the vistasession directory contain mrSESSION.mat
cd(mrSessionDir);

% We need a volume view:
dataType = 'Averages';
setVAnatomyPath(t1file);
hvol = initHiddenGray;

% Set the volume view to the current data type and add the RM model
hvol = viewSet(hvol,'curdt',dataType);

% Load mrVista retinotopy Gray file
rmModel  = dirPth.fmri.vistaGrayFitFile;
hvol     = rmSelect(hvol,1,rmModel);
mmPerVox = viewGet(hvol,'mmpervox');

% White surfaces were used in original code. But since we are using pial
% surface later when exporting parameters from freesurfer to brainstorm
% surfaces, it might be wise to chose pial surface here. Don't think it
% will make much difference because the vertices are same for both white
% and pial surface. Only the coordinate values will change for ex in the
% rois.
hemis = {'lh','rh'};

prfDataFile = dir(fullfile(prfDatamrVistaMatDir,'*.mat'));
load(fullfile(prfDatamrVistaMatDir,prfDataFile.name),'prf_par_exp');

% Loop over all the parameters stored in the exported data file:
prfParamNames = fieldnames(prf_par_exp);


% for n_surf = 1:length(surfaces_to_load)
for h = 1:length(hemis)
    
    % Select hemisphere
    curHemi = hemis{h};
    fprintf('(%s): Loading parameters for %s hemisphere\n',mfilename, curHemi);
    
    % load both pial and white to make midgray
    pialSurfFname = fullfile(fsSurfDir, [curHemi,'.pial']);
    whiteSurfFname = fullfile(fsSurfDir,[curHemi,'.white']);
    
    
    % Load mesh using fs_meshFromSurface, this creates a mrVista compatible
    % mesh. Using 'my own' function that skips the smoothing:
    mrVistaMeshPial = mprf_fs_meshFromSurface(pialSurfFname);
    mrVistaMeshWhite = mprf_fs_meshFromSurface(whiteSurfFname);
    
    mrVistaMeshMidGray_vertices = mean(cat(3, mrVistaMeshPial.vertices, mrVistaMeshWhite.vertices),3);
    
    % compute mapping using mrmMapVerticesToGray (mrmMapGrayToVertices):
    vertex2GrayMap = mrmMapVerticesToGray(mrVistaMeshMidGray_vertices, viewGet(hvol,'nodes'),...
        mmPerVox,3);
    
    % Select the ones to map to mid gray
    verticesToMap = vertex2GrayMap > 0;
    
    % Loop over every prf parameter
    for param = 1:length(prfParamNames)
        
        curParamName = prfParamNames{param};
        if opt.verbose; fprintf('(%s): Exporting %s parameter\n', mfilename, curParamName); end
        
        % preallocate space
        prfFS = nan(size(vertex2GrayMap));
        
        % Select the data that fall within the vertices to map
        curData = prf_par_exp.(curParamName);
        prfFS(verticesToMap) = curData(vertex2GrayMap(verticesToMap));
        
        % store new mapped prf estimates in allData struct
        allData.(curHemi).(curParamName) = prfFS;
        
        % output the results:
        fname = fullfile(prfFSDir,[curHemi '.' curParamName]);
        write_curv(fname,prfFS, 1);
        
    end
    
    
    
end


%% ------------------------------------------------------------------------
% Exporting ROI
%--------------------------------------------------------------------------

% Either export hand-drawn on mrVista mesh to the Freesurfer space,
% which can be set by opt.roimrvToFS = true
% Otherwise locate probabilistic atlas from Wang et al (2015)

if ~opt.roimrvToFS
    fprintf('(%s): Saving a copy of Wang et al. 2015 probabilistic atlas on FreeSurfer surface in roi dir: %s\n', mfilename, roiNames{rr})

    roifname      = @(hem)(sprintf('%s/surf/%s.wang2015_atlas.mgz', dirPth.fs.segPth, hem));
    surfdat       = @(lh, rh, rc)(setfield(setfield([], 'lh', lh.vol(:)), 'rh', rc*rh.vol(:)));
    roiFSWang     = surfdat(MRIread(roifname('lh')), MRIread(roifname('rh')), 1);
    
    for h = 1:length(hemis)
        curHemi = hemis{h};
        outFileName = [curHemi '.wang2015_atlas'];
        fname = fullfile(roiFSDir,outFileName);
        write_curv(fname, roiFSWang.(curHemi),1);
    end
    
    
elseif opt.roimrvToFS
    
    for h = 1:length(hemis)
        curHemi = hemis{h};
        
        if strcmp(curHemi,'lh')
            roisMrV = dir(fullfile(roiDir,'L*.mat'));
            roiNames = {roisMrV.name};
        elseif strcmp(curHemi, 'rh')
            roisMrV = dir(fullfile(roiDir,'R*.mat'));
            roiNames = {roisMrV.name};
        end
        
        
        % Set the roiDilateIterations parameter te prevent 'patchy' outcomes
        prefs = mrmPreferences;
        prefs.roiDilateIterations = 2;
        mrmPreferences(prefs);
        
        % Keep track of the ROI vertices. We de this two times: one for the
        % original ROI vertices, and one for the dilated ROI vertices.
        % Dilation will cause ROIs to claim vertices of adjacent ROIs. We use
        % these variables to keep track of which vertices originally belonged
        % to which ROI and prune the ROIs later on:
        roiVertexIdx_orig = cell(1,length(roiNames));
        roiVertexIdx_dilated = cell(1,length(roiNames));
        
        
        % Loop over the ROIS
        for rr = 1:length(roiNames)
            
            fprintf('(%s): Mapping hand drawn ROI to freesurfer surface: %s\n', mfilename, roiNames{rr})
            
            % Load the current ROI
            load(fullfile(roiDir,roiNames{rr}))
            
            % Get the ROI indices in the gray graph and load the nodes:
            [~, idx] = intersectCols(hvol.coords, ROI.coords);
            roiNodes = hvol.nodes(:,idx);
            
            % load both pial and white to make midgray
            pialSurfFname = fullfile(fsSurfDir, [curHemi,'.pial']);
            whiteSurfFname = fullfile(fsSurfDir,[curHemi,'.white']);
            
            
            % Load mesh using fs_meshFromSurface, this creates a mrVista compatible
            % mesh. Using 'my own' function that skips the smoothing:
            mrVistaMeshPial = mprf_fs_meshFromSurface(pialSurfFname);
            mrVistaMeshWhite = mprf_fs_meshFromSurface(whiteSurfFname);
            mrVistaMeshMidGray_vertices = mean(cat(3, mrVistaMeshPial.vertices, mrVistaMeshWhite.vertices),3);
            
            mrVistaMeshMidGray = mrVistaMeshWhite;
            mrVistaMeshMidGray.vertices = mrVistaMeshMidGray_vertices;
            
            hvol = viewSet(hvol,'add and select mesh',mrVistaMeshMidGray);
            
            
            % Map the gray nodes to the surface vertices, the results are
            % much cleaner than when mapping the other way around.
            vertex2ROIMap = mrmMapGrayToVertices(roiNodes, mrVistaMeshMidGray_vertices, mmPerVox, 2);
            
            % Not all voxels can be mapped, give a warning of more than 25% but
            % less than 95% could not be mapped.
            missing = mean(vertex2ROIMap == 0);
            if missing  > .25 && missing < .95
                warning('Could not map %0.2f percent of voxels to surface for %s', missing *100, roiNames{rr})
            elseif  missing >= .95
                warning('Could not map %0.2f percent of voxels to surface for %s.\nSKIPPING THIS ROI', missing *100, roiNames{rr})
            end
            
            % ROIs can look patchy on a cortical surface, so dilate them a
            % little be to remove holes in the ROI
            % Keep the orignal and dilated ROI indices
            roiVertInds = unique(vertex2ROIMap(vertex2ROIMap > 0));
            roiVertInds2 = adjustPerimeter(roiVertInds, 0, hvol, prefs);
            
            roiVertexIdx_orig{rr} = roiVertInds;
            roiVertexIdx_dilated{rr} = roiVertInds2;
            
        end
        
        fprintf('(%s): Done.\n', mfilename)
        
        % Prune the ROIs to remove overlap due to dilate, keeping the original
        % boundaries between ROIS:
        prunedROIIdx = mprfResolveDilatedROIOverlap(roiVertexIdx_dilated, roiVertexIdx_orig);
        
        % Preallocate space for ROIs mapped to vertices
        allROIs = nan(length(roiNames),length(mrVistaMeshMidGray_vertices));
        
        fprintf('(%s): Combining ROIs for hemi %s \n', mfilename, curHemi)
        
        % Now, store the pruned ROI indices, both as a separate file
        for rr = 1:length(roiNames)
            
            [~,roiName] = fileparts(roiNames{rr});
            
            if ~isempty(prunedROIIdx{rr})
                
                allROIs(rr,prunedROIIdx{rr}) = rr;
                
                outFileName = [curHemi '.' roiName(2:end)];
                fname = fullfile(roiFSDir,outFileName);
                
                write_curv(fname, prunedROIIdx{rr},1);
                
            else
                warning('No indices found for %s, SKIPPING', roiName);
            end
            
        end
        
        fprintf('(%s): Done.\n', mfilename)
        
        % Export a combined ROI
        outFileNameAllROIs = [curHemi '.all_rois'];
        if all(isnan(allROIs(:)))
            warning('The combined ROI could not be created for %s',curHemi)
        else
            fname = fullfile(roiFSDir,outFileNameAllROIs);
            MRIwrite(vol.allROIs,fname);
        end
        
        % Export mask
        mask = false(1,length(mrVistaMeshMidGray_vertices));
        mask(any(allROIs,1))=true;
        outFileNameMask = [curHemi '.mask'];
        
        fname = fullfile(roiFSDir,outFileNameMask);
        write_curv(fname, mask,1);
        
    end
    
end
