%% Configuration

clear('all');
tbUse('vistasoft');

% If these are set, then they won't be deduced from the script location in
% the next section:
sub_name = 'wlsubj030';
session_path = '/Volumes/server/Projects/MEG/Retinotopy/Data/fMRI/wlsubj030/';
session_name = 'MRI_data';
bids_sess    = 'nyu3T01';
stimuli_path = fullfile(session_path, 'Stimuli');
output_path = fullfile(session_path, 'Outputs');

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
    script_path = mfilename('fullpath');
    [code_path,     script_name,  ~] = fileparts(script_path);
    [session_path,  code_name,    ~] = fileparts(code_path);
    [sub_path,      session_name, ~] = fileparts(session_path);
    [subjects_path, sub_name,     ~] = fileparts(sub_path);

    if ~strcmpi(code_name, 'code')
        error(['If initialization script is not in <subj>/<sess>/code' ...
               ' then it must be edited manually to include paths']);
    end
end
if ~exist('session_name', 'var') || ~exist('session_path', 'var')
    error('No session name/path deduced or provided');
end
fprintf('Subject: %-12s  Session: %-20s\n', sub_name, session_name);
   
% Next, figure out freesurfer data if not given
if ~exist('freesurfer_id', 'var'), freesurfer_id = sub_name; end
if ~exist('fs_subjects_dir', 'var')
    fs_subjects_dir = getenv('SUBJECTS_DIR');
end

% Last, figure out the output directory if not given
if ~exist('output_path', 'var')
    output_path = fullfile(session_path, 'Outputs');
end
if ~exist(output_path, 'dir')
    error(sprintf('output_path (%s) not found', output_path));
end


%% Navigate and initialize session

cd(session_path);

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
    functionals2nifti(vw, [], sprintf('%s/%s-xcrds.nii.gz', output_path, outnames{ii}));
    vw = rmLoad(vw, 1, 'y0', 'map');
    functionals2nifti(vw, [], sprintf('%s/%s-ycrds.nii.gz', output_path, outnames{ii}));
    % we also need the variance explained and sigma
    vw = rmLoad(vw, 1, 'sigma', 'map');
    functionals2nifti(vw, [], sprintf('%s/%s-sigma.nii.gz', output_path, outnames{ii}));
    vw = rmLoad(vw, 1, 'variance explained', 'map');
    functionals2nifti(vw, [], sprintf('%s/%s-vexpl.nii.gz', output_path, outnames{ii}));
end

% That's all!