% Ask for input files:
% - Collect the necessary files:
%       1. mrVista directory of subject (directory)
%            In order to run, we need a Volume view and pRF models in the Volume view
%            The anatomy and segmentation for this session have to come
%            directly from a freesurfer session.
%
%       2. mrVista T1. (File)
%            The T1 on which the volume view is based
%
% 		3. mrVista Class file
%            The segmentation of the T1 above
%       
%       4. Brainstorm head model file
%            A head model created in Brainstorm based on coregistration 
%            of subject head with the MEG scanner and the MEG data that we want 
%            to analyse in this mprfSESSION.
%            Importantly, Brainstorm imports a freesurfer directory and the 
%            T1 and segmentation above *have to come from the same freesurfer session*.
%
%       5. Brainstorm anatomy directory
%           Brainstorm imports a freesurfer directory and stores the relevant files
%           and transformed surfaces in a separate directory. We need this one    
%
%       6. Retinotopic model
%           From the mrVista session selected under 1. The pRF models
%           whose parameters we want to use for this session.
%
%       7. Left hemisphere ROIS.
%           I believe the code assumes that left hemisphere ROIs are
%           *prepended* by 'L' or 'l'. Make sure this is the case. As we
%           Brainstorm relies on freesurfer surfaces and it distinguishes
%           between the left and right hemisphere, we *need* the ROIs to be
%           separated by hemisphere as well!
%       8. Right hemisphere ROIS.
%           Same as left hemisphere ROIs, but of course prepended with 'R'
%           or 'r'.
%       
%       9. MEG data
%           A single .sqd file collected for the current subject. You can
%           combine several .sqd files into one, if needed.
%
%       10. MEG stimulus folder
%           A bit cumbersome, but we rely on a stimulus file for every
%           stimulus run (so 19 in total, for both subjects we have
%           thus far). I forgot why. But be sure to select a folder that
%           has all stimulus files, and only the stimulus files, nothing
%           more!
%
%       11. MEG parameter folder
%           While collecting MEG data, we keep track of the stimulus
%           timing. As MEG data is collected on the millisecond scale it is
%           important to check if there are any dropped frames (16 ms
%           delay). So, we need a parameter file for every stimulus run to
%           check for this. We also use these files to match the MEG data
%           to the stimulus presentation. Although the MEG scanner says it
%           has a frequency of 1000 Hz (i.e. 1/0.001), the frequency is actually
%           999.970000899973 Hz (1/0.00100003). In other words a MEG sample
%           is ~31 nanaseconds longer than 1 ms. This does not look like
%           much, but when collecting data for over an hour, this
%           accumulates to a 'delay of the MEG acquisition' of about 138 milliseconds.
%
%       12 MRI stimulus folder
%           Same as for the MEG stimulus folder, but now for MRI.
%           Also forgot why I did it this way.

project_pth = '/Volumes/server/Projects/MEG/Retinotopy/';
subject = 'wlsubj030';
cd(fullfile(project_pth, 'Data', 'fMRI', subject));

global mprfSESSION %#ok<NUSED>

% Check if there is already a mprfSESSION.mat file
if exist(fullfile(pwd,'mprfSESSION.mat'),'file')
    mprf__load_mprfsession('load'); % If so, load it
else
    mprf__load_mprfsession('init'); % If no, create it
    
end

% Initialise the directory structure for the current session:
mprf__init_directory_structure; 

% Start GUI, ask for the location of the files we need, add the paths to
% the mprfSESSION variable:
mprf__session_init_gui;

% We keep track of which paths have been changed, update the paths that
% have been changed and (reimport) the necessary files:
mprf__update_session;

% Import the retinotopic stimulus, i.e. the one that was used to model the
% pRFs in the mrVista session above. It relies on a functioning mrVista
% session.
mprf__import_stimulus('rm_stim');
    
% Import the MEG stimulus, i.e. it will put it in the correct format. It
% will also ask for a grid file (i.e. the meshgrid used to define the
% stimulus), this should actually be part of the stimulus files, but it
% isn't.
mprf__import_stimulus('meg_stim');
    



save('mprfSESSION.mat','mprfSESSION');
clear all
% TODO: make gui, add fields for the new data paths
% Make an update script, laod fields from mprfSESSION.orig.path and put
% them in Gui, check which have changed and update them accordingly.
% Check which fields are available and push data automatically through all
% possible steps. Or allow an optoin like: 'Run','Save session' and thats
% it. I.e. RUn with new fields, or just save the session like that. Then
% make sure analyses can be run on a separate occasion. Make an import
% stimulus GUI -> prediction stimulus, ideally from the imported stimulus
% files and param files. Make sure some info about the stimulus can be
% printed to the Matlab command window, make sure this can be done for both
% the MRI and MEG stimulus. Make a GUI that allows the user to produce a
% predicted MEG signal, make sure this creates a new forlder with all the
% relevant information





















