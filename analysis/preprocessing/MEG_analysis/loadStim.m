function meg_stim = loadStim(subjID, dirPth, opt)
% Function to load stimulus used in the MEG and MRI experiment
%
%   meg_stim = loadStim(subjID, dirPth, opt)
%
% INPUTS:
%   subjID          : subject name (string)
%   dirPth          : paths to files for given subject (struct)
%   opt             : pipeline options (struct with boolean flags).
%
% OUTPUT:
%   meg_stim        : struct with stimulus used in MEG experiment.
%                       Contains X, Y coordinates, rescaled images (for
%                       mrVista solve pRF analysis) and windowed images
%
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019

if opt.verbose; sprintf('(%s) Load MEG stim for subject %s...\n', mfilename, subjID); end

meg_stim = struct();
meg_stim.MEG = loadMEGStimulus(dirPth, opt);
meg_stim.MRI = loadMRIStimulus(dirPth, opt);


end