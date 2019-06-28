function rootPath = mprf_rootPath

% Function to get the current path of this rootpath file. This will be used
% as the base folder of this project.

rootPath=which('mprf_rootPath');

rootPath=fileparts(rootPath);

end






