function mprfSessionInit(do_ernie)
% Use this file to initiliza a mprf session
% It only wants to know if you want to use Ernie's example data or not
% (true/false). Ernie also has brainstorm surfaces and models, needed to
% run a simulation.
%
% prepare your environment by navigating to
% /Volumes/server/Projects/MEG/Retinotopy and add the 'Code' directory to
% your Matlab path:
% addpath(genpath(fullfile(pwd, 'Code')))
%
% This code requires vistasoft to run and will try to add it to your path
% if the vistaRootPath.m file does not exist. It will use tbUse to do so.
% If you do not have tbUse, add vistasoft manually to your path before
% running this script
% If you are using Ernie's data set, you will also need the Remote Data
% Toolbox (RDT).
%
% After initialization, use:
%
% - mprfSessionSmoothExportPRFParams to export and smooth the pRF parameters
%   from the retinotopic model. This file takes one argument.
% 
% - mprfSessionExportPRFDataToFreeSurferSurface to export the pRF data from
%   mrVista format to the freesurfer surface
%
% - mprfSessionExportVistaROIToFreesurferSurface to export ROIs drawn in
%   mrVista to the freesurfer surface
%
% - mprfSessionImportWangAtlas to import the Wang atlas produced by Noah
%   Benson's docker
%
% - mprfSessionFreesurferSurfaceToBrainstorm to transform the freesurface
%   surface prf data and rois to the brainstorm
%
% - mprfSessionRenderDataOnFreesurferSurface to render any of the surface
%  data on a freesurfer surface
%
%
% - mprfSessionRenderDataOnBrainstormSurface to render any of the surface
%   data on a brainstorm surface
%
%

%% Initialize the session

global mprfSESSION
global paths

if exist(fullfile(pwd,'mprfSESSION.mat'),'file')
    answer = questdlg(sprintf('An existing session file has been found. \n Update this file or overwrite'),...
        'Update or overwrite',...
        'Update',...
        'Overwrite',...
        'Update');
    
    if strcmpi(answer,'update')
        load('mprfSESSION.mat');
        mprfSessionFillPathVariable;
        
        
    else

        
    end
    
end



if ~exist('do_ernie','var') || isempty(do_ernie)
    do_ernie = false;
end

if exist('vistaRootPath','file');
    
else
   try
      tbUse('vistasoft');
   catch ME
       error('Vistasoft is not on your Matlab path')
       
   end
    
end

%% update the mprfSESSION variable with the initilization:
mprfUpdateSession(1, 'init');

% Freesurfer directory, if available. Getenv will return empty is variable
% is not found, so that is fine. Could be removed possibly.
paths.orig.fssubjectsdir = getenv('SUBJECTS_DIR');


%% Collect the paths to the data we need. Check for Ernie's data if we are using that:
if do_ernie
    if ~isfield(paths.orig.vista_dir) || isempty(paths.orig.vista_dir)
        % Check if Ernie data is installed:
        paths.orig.vista_dir = mrtInstallSampleData('functional', 'erniePRF',[],false);
    end
    
    if ~isfield(paths.orig.fs_subject_dir) || isempty(paths.orig.fs_subject_dir)
        % Check if anatomy directory is installed:
        paths.orig.fs_subject_dir = mrtInstallSampleData('anatomy/freesurfer', 'ernie', ...
            paths.orig.fssubjectsdir, false);
    end
end

% Root path to the mprfSESSION code directory, similar to vistaRootPath
paths.orig.root_path = mprfRootPath;

% Ask for the remaining paths using a gui:
mprfSessionInit_gui;

% Store the paths:
mprfUpdateSession_Orig_paths;


%% Import the data from its original destiation to the session folder's source data:
mprfSessionImportSourceData;

% Save the mprfSESSION variable:
save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')

clear paths
clear mprfSESSION

end

























