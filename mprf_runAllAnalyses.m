% mprf_runAllAnalysis.m
%
% Master script to run 3 types of analyses of the manuscript: 
%   1. Original: predict MEG responses for single subject, using initally
%      estimated fMRI pRF parameters
%
%   2. Systematically varying pRF size: use predictions from pRFs 
%      where sigma is systematically scaled (up or down) from original 
%      estimate to explain MEG responses.
%
%   3. Systematically varying pRF position: use predictions from pRFs 
%      where [x,y] is systematically rotated around the fovea to explain 
%      MEG responses.

% A population receptive field model of the magnetoencephalography response
% by Kupers, Edadan, Benson, Zuiderbaan, de Jong, Dumoulin and Winawer.
% YEAR. JOURNAL. DOI.
%
% Default script runs the modelfits with 2 free parameters, a gain and a
% reference phase. Model uses split-half cross-validation to get the 
% the final reference phase and betas within a subject. 
%
% To visualize the results, you run the script mprf_makeManuscriptFigures.m
%
% Written by EK @ NYU

%% Define all subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'}; %

% Define additional free params / options here, such as:
saveFig = false;
verbose = false;

%% Fit model with fMRI pRF parameters
opt = getOpts('saveFig', saveFig,'verbose',verbose); % see getOpts function for defaults and more options

for s = length(subjects)
    subjID = subjects{s};
    mprf_main(subjID, opt);
end

%% Fit model with altered pRF size
opt = getOpts('saveFig',saveFig,'verbose',verbose, 'perturbOrigPRFs','size');
for s = 1:length(subjects)
    subjID = subjects{s};
    mprf_main(subjID, opt);
    close all
end

%% Fit model with rotated pRF position
opt = getOpts('saveFig',saveFig,'verbose',verbose, 'perturbOrigPRFs','position');
for s = 1:length(subjects)
    subjID = subjects{s};
    mprf_main(subjID, opt);
    close all
end


%% Fit model with pRF parameters using the group average retinotopy data set collected at NYU 3T scanner.
for s = 1:length(subjects)
    subjID = subjects{s};
    opt = getOpts('saveFig', saveFig,'verbose',verbose, ...
        'useNYU3TAveMaps', true); 
    mprf_main(subjID, opt);
    close all;
end


%% Fit model with scrambled pRF parameters across cortex.
%       opt = getOpts('saveFig',saveFig,'verbose',verbose, 'perturbOrigPRFs','scramble'); 
%       for s = 1:length(subjects)
%           subjID = subjects{s};
%           mprf_main(subjID, opt);
%           close all
%       end


return
