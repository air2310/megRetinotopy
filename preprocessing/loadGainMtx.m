function gainMtx = loadGainMtx(subjID, dirPth, opt)
%
% Function to load gain matrix produced by Brainstorm

if opt.verbose; sprintf('(%s) Load MEG Gain matrix for subject %s...\n', mfilename, subjID); end

if opt.fullSizeGainMtx
    d = dir(fullfile(dirPth.bs.dataPth, '*', 'headmodel*meg_02.mat'));
else
    d = dir(fullfile(dirPth.bs.dataPth, '*', 'headmodel*meg.mat'));
end

gainMtx = load(fullfile(d.folder,d.name));


end