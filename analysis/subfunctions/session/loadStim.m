function meg_stim = loadStim(stimPths)

% Function to load stimulus used in the MEG and MRI experiment

meg_stim = struct();

meg_stim.MEG = loadMEGStimulus(stimPths);
meg_stim.MRI = loadMRIStimulus(stimPths);


end