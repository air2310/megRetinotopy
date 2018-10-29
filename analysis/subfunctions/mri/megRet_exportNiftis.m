function [] = megRet_exportNiftis(subject, sessionDir)

%% Configuration

outputPath = fullfile(sessionDir, 'Outputs');

% If you don't have SUBJECTS_DIR set, then you'll need to set this manually
% to be your FreeSurfer subject's directory
% fs_subjects_dir = '/Volumes/server/Freesurfer_subjects';

% If your subject has a different FreeSurfer subject ID than VistaSoft ID,
% you must set this to the subject's freesurfer id:
% freesurfer_id = 'wl_subj001';

%% Deducing the remaining configuration data
%  We assume that this file is in <something>/<subject>/code/scriptname.m
%  for the purposes of deducing various paths.

% First, get the session path and subject name
if ~exist('sub_name', 'var')
    scriptPath = mfilename('fullpath');
    [codePath,     scriptName,  ~] = fileparts(scriptPath);
    [sessionDir,  codeName,    ~] = fileparts(codePath);
    [subPath,      sessionName, ~] = fileparts(sessionDir);
    [subjectsPath, subject,     ~] = fileparts(subPath);

    if ~strcmpi(codeName, 'code')
        error(['If initialization script is not in <subj>/<sess>/code' ...
               ' then it must be edited manually to include paths']);
    end
end
if ~exist('sessionName', 'var') || ~exist('sessionPath', 'var')
    error('No session name/path deduced or provided');
end
fprintf('Subject: %-12s  Session: %-20s\n', subject, sessionName);
   
% Next, figure out freesurfer data if not given
if ~exist('freesurferID', 'var'), freesurferID = subject; end
if ~exist('fsSubjectsDir', 'var')
    fsSubjectsDir = getenv('SUBJECTS_DIR');
end

% Last, figure out the output directory if not given
if ~exist('outputPath', 'var')
    outputPath = fullfile(sessionDir, 'Outputs');
end
if ~exist(outputPath, 'dir')
    error(sprintf('outputPath (%s) not found', outputPath));
end


%% Navigate and initialize session

cd(sessionDir);

% The datsets we need to deal with:
datasets = {'Averages'};
outnames = {'averages'};

% Do each in turn...
vw = initHiddenGray();
for ii = 1:numel(datasets)
    vw = viewSet(vw, 'current dt', datasets{ii});
    % load a model
    vw = rmSelect(vw, true, sprintf('Gray/%s/rm_%s-fFit.mat', datasets{ii}, datasets{ii}));
    rm = viewGet(vw, 'retinotopy model');

    % Load and export the angle/eccen
    vw = rmLoad(vw, 1, 'x0', 'map');
    vw = viewSet(vw, 'displaymode', 'map');
    functionals2nifti(vw, [], sprintf('%s/%s-xcrds.nii.gz', outputPath, outnames{ii}));
    vw = rmLoad(vw, 1, 'y0', 'map');
    functionals2nifti(vw, [], sprintf('%s/%s-ycrds.nii.gz', outputPath, outnames{ii}));
    % we also need the variance explained and sigma
    vw = rmLoad(vw, 1, 'sigma', 'map');
    functionals2nifti(vw, [], sprintf('%s/%s-sigma.nii.gz', outputPath, outnames{ii}));
    vw = rmLoad(vw, 1, 'variance explained', 'map');
    functionals2nifti(vw, [], sprintf('%s/%s-vexpl.nii.gz', outputPath, outnames{ii}));
end

% That's all!