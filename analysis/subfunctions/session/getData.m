function data = getData(subject)

% Function that loads subject MEG, PRF and Gain matrix data

% NOTE THIS FUNCTION IS STILL UNDER CONSTRUCTIOIN -- DO NOT USE

% data = getData(subject)

% INPUTS
% subject       : subjectID (string)


% check dir
curDir = pwd;
[path, subjectName] = fileparts(curDir);
if ~strcmp(subjectName, subject)
    cd(fullfile(mprfRootPath, 'data', 'subjectSession', subject))
end

%% Get MEG data

% Go to dir
dataDir = fullfile('data', 'meg', 'preproc', 'pp');

% Get files
d1 = dir(fullfile(dataDir, '*data*.mat'));
d2 = dir(fullfile(dataDir, '*conditions*.mat'));

% Load file
data.meg = load(fullfile(d1.folder,d1.name));
data.meg_conditions = load(fullfile(d2.folder,d2.name));

%% Get prf data

% go to dir

% get file

% load file

%% Get gain matrix data
% go to dir

% get file

% load file