function roiData = loadWangROIs(roiDir, searchName, combineHemi, combineRois)

% Find ROI file from roiDir with searchName
roiFile = dir(fullfile(roiDir,searchName));

% Define ROI names
wangROINames = {'V1v', 'V1d', 'V2v', 'V2d', 'V3v', 'V3d', 'hV4', 'VO1', 'VO2', 'PHC1', 'PHC2', ...
    'TO2', 'TO1', 'LO2', 'LO1', 'V3B', 'V3A', 'IPS0', 'IPS1', 'IPS2', 'IPS3', 'IPS4', ...
    'IPS5', 'SPL1', 'FEF'};

if combineRois
    roisToCombine = {'V1','V2','V3','IPS','LO','TO','VO', 'PHC'};
    roisIndiv     = find(~contains(wangROINames,roisToCombine));
end


% Create struct for roi data output
roiData = struct();

% Loop over the LH roi files
curFile = roiFile.name;
tmp = str_split(curFile, '.');

if strcmp(tmp{2},'mat')
    error('(%s): This function does not support mrVista matfiles yet', mfilename)
else
    roiName = tmp{2};
end

% Find the corresponding rh file:
if combineHemi
    if strcmp(tmp{1}, 'lh')
        h = 'rh.';
    elseif strcmp(tmp{1}, 'rh')
        h = 'lh.';
    else
        error('(%s): current roi %s does not seem to have other hemi file', mfilename, tmp{2})
    end
    
    if strcmp(tmp{3}, 'mgz')
        otherHemiFile = [h roiName '.' tmp{3}];
    else
        otherHemiFile = strcat(h,roiName);
    end
                
    % load and concatenate:
    vert1 = MRIread(fullfile(roiDir,curFile));
    vert2 = MRIread(fullfile(roiDir,otherHemiFile));
    roiVertices = [squeeze(vert1.vol); squeeze(vert2.vol)];
else
    vert = read_curv(fullfile(roiDir,curFile));
    roiVertices = vert;
end

% Add all ROIs and mask
roiData.allROIs = roiVertices;
roiData.mask    = roiVertices>0;
roiData.V123    = zeros(size(roiVertices));
roiData.V123(ismember(roiVertices,1:6)) = find(ismember(roiVertices,1:6));
roiData.V123mask = false(size(roiVertices));
roiData.V123mask(find(roiData.V123)) = true;

% Combine data of certain ROIs if requested
if combineRois
    
    for rc = 1:length(roisToCombine)
        match = find(contains(wangROINames,roisToCombine(rc)));
        
        tmpData = [];
        
        for ii = 1:length(match)
            indices = find(roiVertices==match(ii));
            tmpData = [tmpData; indices];
        end
        
        roi = zeros(size(roiVertices));
        tmpIdx = unique(tmpData, 'rows');
        roi(tmpIdx)=tmpIdx;
        roiData.(roisToCombine{rc}) = roi;
    end
    
    for rc = roisIndiv
        roi = zeros(size(roiVertices));
        roi(roiVertices==rc) = rc;
        roiData.(wangROINames{rc}) =  roi;
    end
    
else % otherwise just add 
    
    for ii = 1:length(wangROINames)
        % Add ROIs that do not need to be combined
        roi = zeros(size(roiVertices));
        roi(roiVertices==ii) = ii;
        roiData.(wangROINames{ii}) = roi;
    end
end

return