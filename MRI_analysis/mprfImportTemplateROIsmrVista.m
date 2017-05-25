% Run from subject directory

mrVista('3');


VOLUME{1} = meshBuild(VOLUME{1}, 'left');
VOLUME{1} = meshBuild(VOLUME{1}, 'right');


% Inflate the meshes
msh = cell(1,2);
for ii = 1:2
    msh{ii} = viewGet(VOLUME{1}, 'Mesh', ii);
    msh{ii} = meshSet(msh{ii}, 'smooth iterations',200);
    msh{ii} = meshSet(msh{ii}, 'smooth relaxation',1);
    msh{ii} = meshSmooth(msh{ii});
end
VOLUME{1} = viewSet(VOLUME{1}, 'all meshes', msh);


% Select wang atlas niftis
[fname, pth] = uigetfile('*.mgz','Select wang atlas file');




cur_path = fullfile(pth, fname);


ni = MRIread(cur_path);

nifti_path = fullfile(pth, sprintf('%s.nii.gz', fname));
ni.vol = abs(ni.vol);
MRIwrite(ni, nifti_path);


hideGrayROIs = viewGet(VOLUME{1}, 'hide gray ROIs');
VOLUME{1} = viewSet(VOLUME{1}, 'hide Gray ROIs', true);


Wang_ROI_Names = {'V1' 'V2' 'V3' 'hV4', 'VO1' 'LO1' 'V3A' 'V3B'};

wangComments = 'ROI made based on maximum probabilty map from Wang et al, 2014 ''Probabilistic Maps of Visual Topography in Human Cortex''';


% count the current number of ROIs
nROIs = length(viewGet(VOLUME{1}, 'ROIs'));

% create new ROIs from atlas
VOLUME{1} = nifti2ROI(VOLUME{1}, nifti_path);

% how many new ROIs did we get?
numNewROIs = length(viewGet(VOLUME{1}, 'ROIs')) - nROIs;
assert(numNewROIs>0)

for jj = 1:numNewROIs
    this_roi       = viewGet(VOLUME{1}, 'ROI name', nROIs + jj);
    this_roi_num   = str2double(this_roi(regexp(this_roi, '[0-9]')));
    this_roi_name  = sprintf('WangAtlas_%s', Wang_ROI_Names{this_roi_num});
    VOLUME{1}             = viewSet(VOLUME{1}, 'ROI Name', this_roi_name, nROIs + jj);
    comments       = viewGet(VOLUME{1}, 'ROI comments', nROIs + jj);
    comments       = sprintf('%s\n%s', comments, wangComments);
    VOLUME{1}             = viewSet(VOLUME{1}, 'ROI comments', comments, nROIs + jj);
end



VOLUME{1} = viewSet(VOLUME{1}, 'hide Gray ROIs', hideGrayROIs);








% Save the ROIs
local = true; forceSave = true;
saveAllROIs(VOLUME{1}, local, forceSave);

% Let's look at the ROIs on meshes
%   Store the coords to vertex mapping for each ROI for quicker drawing
VOLUME{1} = roiSetVertIndsAllMeshes(VOLUME{1});

VOLUME{1} = meshUpdateAll(VOLUME{1});






















