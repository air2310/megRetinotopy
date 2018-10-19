function [] = initVista(subject, scriptPth, sessionFolder, bidsSession)

%% How this file works:
%  (1) You create a VistaSoft subject directory (but not an anatomy
%      directory)
%  (2) You process your subject with FreeSurfer
%  (3) You make Raw, Preproc, and Code directories in your VistaSoft dir
%  (4) You download your raw data to the Raw directory and run it with
%      Serra's preproc python script (such that the Preproc directory is
%      populated)
%  (5) You run this script to initialize anatomy and VistaSoft data

%% Configuration

% This directory should hold the path to the anatomies for individual
% subjects:
anatomyPath = '/Volumes/server/Projects/Anatomy';

% If you don't have SUBJECTS_DIR set, then you'll need to set this manually
% to be your FreeSurfer subject's directory
freesurferSubjectDir = '/Volumes/server/Freesurfer_subjects';

% If your subject has a different FreeSurfer subject ID than VistaSoft ID,
% you must set this to the subject's freesurfer id:
freesurferID = subject;

% The annot_pattern is the string used to create the scan names via the
% sprintf function; e.g., if you have 4 scans that are all pRF scans, and
% you choose 'myPRF%02d' for this string,
% then your scans will be named myPRF01, myPRF02, myPRF03, and myPRF04.
annot_pattern = 'PRF%02d';

%% Deducing the remaining configuration data
%  We assume that this file is in <something>/<subject>/code/scriptname.m
%  for the purposes of deducing various paths.

% First, get the session path and subject name
if ~exist('subject', 'var')
    scriptPath = mfilename('fullpath');
    disp(scriptPath);
    [codePath,     scriptName,  ~] = fileparts(scriptPath);
    [sessionPth,  codeName,    ~]  = fileparts(codePath);
    [subPath,      sessionName, ~] = fileparts(sessionPth);
    [subjectsPath, subject,     ~] = fileparts(subPath);
    
    if ~strcmpi(codeName, 'code')
        error(['If initialization script is not in <subj>/<sess>/code'
            ' then it must be edited manually to include paths']);
    end
end
if ~exist('sessionPame', 'var') || ~exist('sessionPath', 'var')
    error('No session name/path deduced or provided');
end
fprintf('Subject: %-12s  Session: %-20s\n', subject, sessionName);

% Next, figure out freesurfer data if not given
if ~exist('freesurferID', 'var'), freesurferID = subject; end
if ~exist('freesurferSubjectDir', 'var')
    freesurferSubjectDir = getenv('SUBJECTS_DIR');
end

%% Initializing the 3DAnatomy directory

% Whole brain reference
anatomyPath = fullfile(sessionPath, '3DAnatomy');
if exist(anatomyPath, 'dir')
    fprintf('Found anatomy directory %s...\n', anatomyPath);
else
    extAnatomyPath = fullfile(anatomyPath, subject);
    if ~exist(extAnatomyPath, 'dir')
        fprintf('Initializing anatomy data in %s...\n', extAnatomyPath);
        init_anatomy(extAnatomyPath, freesurferID, freesurferSubjectDir);
    end
    fprintf('Linking anatomy directory %s to %s...\n', ...
        extAnatomyPath, anatomy_path);
    cmd = sprintf('ln -s %s %s', extAnatomyPath, anatomy_path);
    [sysStatus, sysRes] = system(cmd);
    if sysStatus ~= 0
        error(sprintf('Could not link anatomy dir %s to %s;\n%s\n', ...
            extAnatomyPath, anatomyPath, sysRes));
    end
end

t1File = fullfile(anatomyPath, 't1.nii.gz');
if ~exist(t1File,'file'), error('Anatomy file %s not found', t1File); end


%% Finding the EPI files

% Find the preproc path...
preprocPaths = {'preproc', 'Preproc', 'preprocessed', 'Preprocessed'};
for i = 1:numel(preprocPaths)
    preprocPath = fullfile(sessionPath, sessionName, 'raw', 'derivatives', preprocPaths{i});
    if isdir(preprocPath), break;
    else                    preprocPath = [];
    end
end
if isempty(preprocPath), error('Could not find a preproc directory');
else fprintf('Using preproc directory %s...\n', preprocPath);
end

d = dir(fullfile(preprocPath, sprintf('sub-%s', subject), sprintf('ses-%s',bidsSession), '*preproc.nii.gz'));
for ii = 1:length(d)
    epiFiles{ii} = fullfile(preprocPath, sprintf('sub-%s', subject), sprintf('ses-%s',bidsSession), d(ii).name);
end
fprintf('Found %d preprocessed EPI files...\n', numel(epiFiles));

% Inplane (mean of unwarped, merged, distortion scan)
ip = fullfile(preprocPath, sprintf('sub-%s', subject), sprintf('ses-%s',bidsSession), 'distortion_merged_corrected_mean.nii.gz');
if ~exist(ip, 'file'), error('Inplane file %s not found', ip); end

% alignment file
alignFile = fullfile(preprocPath, sprintf('sub-%s', subject), sprintf('ses-%s',bidsSession),  'distort2anat_tkreg.dat');
if ~exist(alignFile, 'file')
    error(sprintf('Alignment file %s not found', alignFile));
end



%% Generating the params structure

% trim off trailing / of the session-path if need-be:
while sessionPath(end) == filesep, sessionPath = sessionPath(1:end-1); end

% move ourselves over to the session's directory
cd(sessionPath);

% remove absolute paths where appropriate so that the pRF solutions can be run
% on the HPC
spl = numel(sessionPath);
if startsWith(ip, sessionPath)
    ip = ip(spl+2:end);
end
for ii = 1:numel(epiFiles)
    if startsWith(epiFiles{ii}, session_path)
        epiFiles{ii} = epiFiles{ii}(spl+2:end);
    end
end
if startsWith(t1File, session_path)
    t1File = t1File(spl+2:end);
end

disp(ip);
disp(t1File);

params = mrInitDefaultParams;

% And insert the required parameters:
params.inplane      = ip;
params.functionals  = epiFiles;
params.sessionDir   = session_path;

% Specify some optional parameters
params.vAnatomy = t1File;
params.subject  = subject;

% Name each scan
params.annotations = cell(numel(epiFiles), 1);
for ii = 1:numel(epiFiles)
    params.annotations{ii} = sprintf(annot_pattern, ii);
end

%% Running the initialization
ok = mrInit(params);

% open a session without GUI
userpathvw = initHiddenInplane();

fid = fopen(alignFile, 'r');

% skip header lines
for k=1:4, tline = fgets(fid); end
for k = 1:4; R(k,:) = str2num(fgetl(fid)); end
fclose(fid);
vistaAlignment = fs_tkreg2vista(R, ip, t1File);
mrSESSION = sessionSet(mrSESSION, 'alignment', vistaAlignment);
saveSession();

% Install gray segmentation
% path to class file
t1classPath = fullfile(anatomyPath, 't1_class.nii.gz');
% number of layers in gray graph along surface (3 layers for 1 mm voxels)
numGrayLayers = 3;
installSegmentation([], [], t1classPath, numGrayLayers);

% Import freesurfer meshes
lsrf_file = fullfile(anatomyPath, 'Left', '3DMeshes', 'Left_white.mat');
if ~exist(lsrf_file, 'file')
    fprintf('Initializing FreeSurfer meshes...\n');
    subFSPath = fullfile(freesurferSubjectDir, freesurferID);
    meshImportFreesurferSurfaces(subFSPath, 'b');
end

%% Helper Functions

% Helper function for initializing anatomy from FreeSurfer
    function init_anatomy(anatPath, freesurferID, freesurferSubjectDir)
        % see if we have a subjects dir
        if isempty(freesurferSubjectDir)
            error('no SUBJECTS_DIR given');
        elseif ~exist(freesurferSubjectDir, 'dir')
            error(sprintf('SUBJECTS_DIR not found: %s', freesurferSubjectDir));
        end
        dataDir = fullfile(freesurferSubjectDir, freesurferID);
        if ~exist(dataDir, 'dir')
            error(sprintf('FreeSurfer directory (%s) not found', dataDir));
        end
        % Make the VistaSoft anatomy directory and check it
        mkdir(anatPath);
        if ~exist(anatPath, 'dir')
            error(sprintf('could not make anatomy directory %s', anatPath));
        end
        outfile = fullfile(anatPath, 't1_class.nii.gz');
        t1file = fullfile(anatPath,  't1.nii.gz');
        alignTo = fullfile(dataDir,  'mri', 'orig.mgz');
        fs_ribbon2itk(freesurferID, outfile, true, alignTo);
        % Check that you created a t1 class file (ribbon) and t1 anatomy
        if ~exist(outfile, 'file'),
            error(sprintf('Failed to create class file %s', outfile));
        end
        if ~exist(t1file, 'file'),
            error(sprintf('Failed to create T1 file %s', t1file));
        end
    end

end