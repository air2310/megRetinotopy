% Ask for input files:
% - Collect the necessary files:
% 		1. Retinotopic model
% 			- Also import retinotopic stimulus
%           - Keep this stimulus separate from other stimuli
%           - Import stimulus using a GUI, allow for annotation makes
%           selecting the correct stimulus easier. So create a GUI to
%           select a stimulus that displays the stimulus settings etc
% 		2. mrVista T1
% 		3. mrVista Class file
% 		4. Brainstorm anatomy folder
% 		5. Brainstorm head model file
% 		6. subject Freesurfer directory
% 		7. mrVista LHS ROIs
% 		8. mrVista RHS ROIs
% 		9. MEG stimulus
% 			- Change stimulus creation logic, i.e. add grid file to the stimulus file -> not now, but perhaps later
% 			- Just import the Stimulus folder in the subjects data directory, assume the following structure
% 				- Stimulus_files
% 				- Param_files
% 				- Other_files
% 				
% 			Or better: create it on initialization
% 		
% 		10. Param file directory (needed for preprocessing the MEG data)
% 		11. Stimulus file directory (needed for preprocessing the MEG data)
% 			- They are not that big, so we can just copy them into the subject folder
% 		
% 		Make sure any field can be left empty. The associated processes just wont run then.
% 		User should be able to pick up and continue processing at any point
% 		
% 		
% 	- Initialize folder structure:
% 		- source
% 		- rois
% 		- prf_data
% 		- predictions
%       - data
% 		
% 		with all the relevant sub_folders: check existing folder
% 		
% 	- Keep track of where to put which files in a separate function, sort of a linker function, similar to mprfExportPath, but then just a little less tedious (hopefully)
% 		
% 
% 

global mprfSESSION %#ok<NUSED>

if exist(fullfile(pwd,'mprfSESSION.mat'),'file')
    mprf__load_mprfsession('load');
else
    mprf__load_mprfsession('init');
    
end
% make sure mprfSESSION exists is by the time the code gets here.
% Make sure it is a global before calling this function, return the paths
% that are created.

mprf__init_directory_structure;
mprf__session_init_gui;

mprf__update_session;


mprf__import_stimulus('rm_stim');
    
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





















