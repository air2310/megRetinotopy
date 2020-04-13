function ok = mprf_createSymbolicLinksToData(pathsToLink)

% Check input argument, if nonexisting, create default
if ~exist('pathsToLink','var') || isempty(pathsToLink, 'var')
    pathsToLink.data = '/Volumes/server/Projects/MEG/Retinotopy';
    pathsToLink.Freesurfer = '/Volumes/server/Freesurfer_Subjects';
    pathsToLink.brainstorm = '/Volumes/server/Projects/MEG/brainstorm_db';
end

% Create data dir in git repo (this folder will be ignored by git)
if ~exist(fullfile(mprf_rootPath,'data'),'dir')
    mkdir(fullfile(mprf_rootPath,'data'))
end


dataDir = filestrp(pathsToLink.data);
FreesurferDir = filestrp(pathsToLink.data);
brainstormDir = filestrp(pathsToLink.brainstorm);

% Check if folders exist, if not, go to data folder to make symbolic link
if ~exist(fullfile(mprf_rootPath,'data', dataDir), 'dir')
    cd(fullfile(mprf_rootPath,'data'))
    system(sprintf('ln -s %s',  pathsToLink.data));
end

if ~exist(fullfile(mprf_rootPath,'data', FreesurferDir), 'dir')
    % Go to data folder to make symbolic links
    cd(fullfile(mprf_rootPath,'data'))
    system(sprintf('ln -s %s',  pathsToLink.Freesurfer));
end

if ~exist(fullfile(mprf_rootPath,'data', brainstormDir), 'dir')
    % Go to data folder to make symbolic links
    cd(fullfile(mprf_rootPath,'data'))
    system(sprintf('ln -s %s',  pathsToLink.brainstorm));
end

return