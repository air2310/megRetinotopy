function [] = megRet_solvePRFs(subject, sessionDir, sessionName)

% If you don't have SUBJECTS_DIR set, then you'll need to set this manually
% to be your FreeSurfer subject's directory
fsSubjectsDir = '/Volumes/server/Freesurfer_subjects';

% If your subject has a different FreeSurfer subject ID than VistaSoft ID,
% you must set this to the subject's freesurfer id:
freesurferID = subject;

% The annot_pattern is the string used to create the scan names via the
% sprintf function; e.g., if you have 4 scans that are all pRF scans, and
% you choose 'myPRF%02d' for this string,
% then your scans will be named myPRF01, myPRF02, myPRF03, and myPRF04.
annotPattern = 'PRF%02d';

%% Deducing the remaining configuration data
%  We assume that this file is in <something>/<subject>/code/scriptname.m
%  for the purposes of deducing various paths.

% First, get the session path and subject name
if ~exist('subject', 'var')
    script_path = mfilename('fullpath');
    [codePath,     scriptName,  ~] = fileparts(script_path);
    [sessionDir,    codeName,    ~] = fileparts(codePath);
    [subPath,       sessionName, ~] = fileparts(sessionDir);
    [subjectsPath,  subject,     ~] = fileparts(subPath);

    if ~strcmpi(codeName, 'code')
        error(['If initialization script is not in <subj>/<sess>/code' ...
               ' then it must be edited manually to include paths']);
    end
end
if ~exist('sessionName', 'var') || ~exist('sessionDir', 'var')
    error('No session name/path deduced or provided');
end
fprintf('Subject: %-12s  Session: %-20s\n', subject, sessionName);
   
% Next, figure out freesurfer data if not given
if ~exist('freesurferID', 'var'), freesurferID = subject; end
if ~exist('fs_subjectsDir', 'var')
    fsSubjectsDir = getenv('SUBJECTS_DIR');
end

% Figure out the stimulus path if not given
if ~exist('paramsFile', 'var') || ~exist('imagesFile', 'var')
    needStimdir = true;
else
    needStimdir = false;
end
if ~exist('stimuliPath', 'dir') && needStimdir
    stimuliPaths = {'Stimulus', 'stimulus', 'Stimuli', 'stimuli', ...
                     'Stim', 'stim'};
    for i = 1:numel(stimuliPaths)
        stimuliPath = fullfile(sessionDir, stimuliPaths{i});
        if isdir(stimuliPath), break;
        else stimuliPath = [];
        end
    end
    if isempty(stimuliPath)
        error('Could not deduce the stimuli path')
    end
end
% And figure out the actual param and image files...
if ~exist('paramsFile', 'var')
    fls = dir(fullfile(stimuliPath, '*_params.mat'));
    paramsFile = [];
    for ii = 1:numel(fls)
        fldat = fls(ii);
        fl = fullfile(stimuliPath, fldat.name);
        if fl(1) == '.',  continue;
        elseif isdir(fl), continue;
        else              paramsFile = fl;
                          break;
        end
    end
    if isempty(paramsFile)
        error('Could not deduce location of parameters mat file');
    end
end
if ~exist(paramsFile, 'file')
    error(sprintf('params file (%s) does not exist', paramsFile));
end
if ~exist('imagesFile', 'var')
    fls = dir(fullfile(stimuliPath, '*_images.mat'));
    imagesFile = [];
    for ii = 1:numel(fls)
        fldat = fls(ii);
        fl = fullfile(stimuliPath, fldat.name);
        if fl(1) == '.',  continue;
        elseif isdir(fl), continue;
        else              imagesFile = fl;
                          break;
        end
    end
    if isempty(imagesFile)
        error('Could not deduce location of images mat file');
    end
end
if ~exist(imagesFile, 'file')
    error(sprintf('images file (%s) does not exist', imagesFile));
end


%% Navigate

cd(sessionDir);
load mrSESSION
gr = initHiddenGray();
gr = viewSet(gr, 'current dt', 'Averages');


%% Initialize Retinotopy Parameters
prfModels = 'one gaussian';
% load(params_file, 'params');

% Set default retinotopy stimulus model parameters
sParams = rmCreateStim(gr);
sParams.stimType   = 'StimFromScan'; % This means the stimulus images will
                                     % be read from a file.
if strcmp(subject, 'wlsubj030')
    sParams.stimSize   = 11.1859;        % usually 12.4; but for experiment the stimsize is limited by MEG screen fov           % stimulus radius (deg visual angle)
elseif strcmp(subject, 'wlsubj068')
    sParams.stimSize   = 10;
end

sParams.nDCT       = 1;              % detrending frequeny maximum (cycles
                                     % per scan): 1 means 3 detrending
                                     % terms, DC (0 cps), 0.5, and 1 cps
sParams.imFile     = imagesFile;    % file containing stimulus images
sParams.paramsFile = paramsFile;    % file containing stimulus parameters
sParams.nCycles     = 1;
sParams.framePeriod = 1.5;
sParams.prescanDuration = 0;


% 'thresholdedBinary',  whenreading in images, treat any pixel value
% different from background as a 1, else 0
sParams.imFilter   = 'thresholdedBinary';
% we switch from the default positive Boynton hRF to the biphasic SPM style
sParams.hrfType    = 'two gammas (SPM style)';
% pre-scan duration will be stored in frames for the rm, but was stored in
% seconds in the stimulus file
sParams.prescanDuration = sParams.prescanDuration/sParams.framePeriod;
saveSession();


%% Solve the retinotopy models
dt = viewGet(gr, 'current dt');
dataTYPES(dt).retinotopyModelParams = [];
dataTYPES(dt) = dtSet(dataTYPES(dt), 'rm stim params', sParams);
saveSession();

gr = rmMain(gr, [], 'coarse to fine', ...
            'model', {prfModels}, ...
            'matFileName', sprintf('rm_%s', dataTYPES(dt).name));


