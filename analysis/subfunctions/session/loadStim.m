function retStim = loadStim(stimPths)

% Function to load stimulus used in the MEG and MRI experiment

retStim = struct();

retStim.MEG = loadMEGStimulus(stimPths);
retStim.MRI = loadMRIStimulus(stimPths);


end