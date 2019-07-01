function meg_stim = loadStim(subjID, dirPth, opt)
%
% Function to load stimulus used in the MEG and MRI experiment

if opt.verbose; sprintf('(%s) Load MEG stim for subject %s...\n', mfilename, subjID); end

meg_stim = struct();
meg_stim.MEG = loadMEGStimulus(dirPth, opt);
% meg_stim.MRI = loadMRIStimulus(dirPth);


end