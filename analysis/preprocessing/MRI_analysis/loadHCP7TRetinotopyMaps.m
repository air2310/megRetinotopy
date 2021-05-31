function prfDataPath = loadHCP7TRetinotopyMaps(subjID, dirPth, opt)
% Function to get pRF data from HCP 7T retinotopy dataset 
% analyzed by Benson et al. (2018), see osf.io/bw9ec/wiki/home
% This script relies on steps subject having average HCP subject data 
% (999999) interpolated to subject's FreeSurfer surface. See hcpave_to_sub.py 
%
%   [prfData, prfDataPath] = loadHCP7TRetinotopyMaps(subjID, dirPth, opt)
%
% INPUTS:
%   subjID          : subject name (string)
%   dirPth          : paths to files for given subject (struct)
%   opt             : pipeline options (struct with boolean flags).
%
% OUTPUTS:
%   prfData         : prf data from average HCP retinotopy interpolated 
%                       back to subject's surface (struct)
%   prfDataPath     : path where prfData are saved (string)
%
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2021



% Check if freesurfer matlab toolbox is added
if ~exist('MRIread', 'file')
    fshome = getenv('FREESURFER_HOME');
    if isempty(fshome)
        setenv('FREESURFER_HOME', '/Applications/freesurfer');
        fshome = getenv('FREESURFER_HOME');
    end
    addpath(genpath(fullfile(fshome,'matlab')));
    addpath(genpath(fullfile(fshome,'fsfast','toolbox')));
end


% create folder to save figures
saveDir  = fullfile(dirPth.fmri.saveDataPth,'HCP7TAveRetinotopyPred');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end
if ~exist(fullfile(saveDir,'prfFS'),'dir'); mkdir(fullfile(saveDir,'prfFS')); end
if ~exist(fullfile(saveDir,'prfBS'),'dir'); mkdir(fullfile(saveDir,'prfBS')); end


if opt.fullSizeMesh
    prfDataPath = fullfile(saveDir,'prfFS');
else
    prfDataPath = fullfile(saveDir,'prfBS');
end

if ~exist(fullfile(prfDataPath,'pial.x'))

    %% 2. Load retinotopy templates
    origPth = fullfile(dirPth.fmri.dataPth, subjID, 'hcpave_interp');

    surfdataFS = struct();

    hemis = {'lh','rh'};
    prfParamNames = {'varexplained', 'beta', 'x', 'y', 'sigma', 'eccentricity','polar_angle'};

    for h = 1:length(hemis)

        % Select hemisphere
        curHemi = hemis{h};
        fprintf('(%s): Loading parameters for %s hemisphere\n',mfilename, curHemi);

        % Loop over every prf parameter
        for param = 1:length(prfParamNames)

            curParamName = prfParamNames{param};
            if opt.verbose; fprintf('(%s): Exporting %s parameter\n', mfilename, curParamName); end

            % Select the data that fall within the vertices to map
            tmp = MRIread(fullfile(origPth,[curHemi '.' curParamName '_output-file.mgz']));
            prfFS = squeeze(tmp.vol);
            
            if strcmp(curParamName,'varexplained')
                prfFS = prfFS./100; % convert from percent to fraction
            end
            
            if strcmp(curParamName,'polar_angle')
                prfFS = deg2rad(prfFS); % convert from degrees to radians
            end
            
            % store new mapped prf estimates in allData struct
            surfdataFS.(curHemi).(curParamName) = prfFS;

            % output the results:
            fname = fullfile(saveDir,'prfFS',[curHemi '.' curParamName]);
            write_curv(fname,prfFS, 1);      
        end  
    end

    % Open up brainstorm if database can't be found
    if ~exist('GlobalData.DataBase.iProtocol', 'var')
        brainstorm nogui;
    end

    % Now combine data
    for param = 1:length(prfParamNames)
        curParamName = prfParamNames{param};

        % For freesurfer
        bothHemiFSData = [surfdataFS.lh.(curParamName); surfdataFS.rh.(curParamName)];

        % Find corresponding BS vertices for FS surface data
        bothHemiBSData = tess_fs2bst(subjID, dirPth.fs.segPth, surfdataFS.lh.(curParamName), surfdataFS.rh.(curParamName));

        % Save the BS results:
        curFileToSave = ['pial.' curParamName];
        fname = fullfile(saveDir,'prfBS',curFileToSave);
        write_curv(fname,bothHemiBSData,1);
        fprintf('(%s): Brainstorm combined hemi files: %s\n',mfilename, curFileToSave);

        % Save the FS results:
        curFileToSave = ['pial.' curParamName];
        fname = fullfile(saveDir,'prfFS',curFileToSave);
        write_curv(fname,bothHemiFSData,1);
        fprintf('(%s): Freesurfer combined hemi files: %s\n',mfilename, curFileToSave);

    end
end

if ~exist(fullfile(prfDataPath,'pial.mask'))
    
    if ~exist(fullfile(saveDir,'roiBS'), 'dir') || ~exist(fullfile(saveDir,'roiFS'), 'dir')
        mkdir(fullfile(saveDir,'roiBS'));
        mkdir(fullfile(saveDir,'roiFS'));
    end
    
    surfaceToMapTo    = 'pial';
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
    write_curv(fullfile(saveDir,'roiBS',curFileToSave),bothHemiBSROI,1);
    write_curv(fullfile(saveDir,'prfBS',curFileToSave),bothHemiBSROI,1);
    fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, curFileToSave);
    
    % COMBINE HEMI FREESURFER file
    write_curv(fullfile(saveDir,'roiFS',curFileToSave),bothHemiFSROI,1);
    write_curv(fullfile(saveDir,'prfFS',curFileToSave),bothHemiFSROI,1);
    fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, curFileToSave);
    
    
    % ---- Store mask of Wang atlas ----
    % BRAINSTORM Pial file
    mask = bothHemiBSROI>0;
    curFileToSave = [surfaceToMapTo '.mask'];
    write_curv(fullfile(saveDir,'roiBS',curFileToSave),mask,1);
    write_curv(fullfile(saveDir,'prfBS',curFileToSave),mask,1);
    fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, curFileToSave);
    
    % COMBINE HEMI FREESURFER file
    mask = bothHemiFSROI>0;
    curFileToSave = [surfaceToMapTo '.mask'];
    write_curv(fullfile(saveDir,'roiFS',curFileToSave),mask,1);
    write_curv(fullfile(saveDir,'prfFS',curFileToSave),mask,1);
    fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, curFileToSave);
    
    % ---- Store just an roi mask of V1-V3 areas in Wang atlas ----
    V123idx = 1:6; % first six rois are "V1v" "V1d" "V2v" "V2d" "V3v" "V3d"
    
    % BRAINSTORM Pial file
    BS_V123mask = ismember(bothHemiBSROI, V123idx);
    curFileToSave = [surfaceToMapTo '.V123mask'];
    write_curv(fullfile(saveDir,'roiBS',curFileToSave),BS_V123mask,1);
    write_curv(fullfile(saveDir,'prfBS',curFileToSave),BS_V123mask,1);
    fprintf('(%s): Brainstorm combined roi files: %s\n',mfilename, curFileToSave);
    
    % COMBINE HEMI FREESURFER file
    FS_V123mask = ismember(bothHemiFSROI, V123idx);
    write_curv(fullfile(saveDir,'roiFS',curFileToSave),FS_V123mask,1);
    write_curv(fullfile(saveDir,'prfFS',curFileToSave),FS_V123mask,1);
    fprintf('(%s): Freesurfer combined roi files: %s\n',mfilename, curFileToSave);
    
end
return
