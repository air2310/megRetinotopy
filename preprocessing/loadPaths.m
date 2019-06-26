function dirPth = loadPaths(subjID)

dirPth  = struct();

dirPth.subjID = subjID;

%% ------ MEG -------
dirPth.meg.dataPth   = './Data/MEG/';
dirPth.meg.saveFigPth    = fullfile('./Quality_check/',subjID, 'meg');

% Derive other file paths
dirPth.meg.rawSqdPath = fullfile(dirPth.meg.dataPth, subjID, 'raw');
dirPth.meg.paramFilePth = fullfile(dirPth.meg.dataPth, subjID,  'paramFiles');
dirPth.meg.stimFilePth = fullfile(dirPth.meg.dataPth, subjID, 'stimFiles');



