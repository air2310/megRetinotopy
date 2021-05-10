% mprf_runAllAnalysis.m
%
% Master script to run 3 types of analyses of the manuscript:
%
% A population receptive field model of the magnetoencephalography response
% by Kupers, Edadan, Benson, Zuiderbaan, de Jong, Dumoulin and Winawer.
% YEAR. JOURNAL. DOI.
%
% Default script runs the modelfits with 2 free parameters, a gain and a
% reference phase. Model uses split-half cross-validation to get the 
% the final reference phase and betas within a subject.
%
% Previous versions of the model could add an additional offset, 
% ('addOffsetParam') or refits gain parameters (so only cross-validate
% reference phase, 'refitGainParam'). But this alternative model fit can be
% subject to overfitting and is therefore not recommended. 
%
% To visualize the results, you run the script mprf_makeManuscriptFigures.m
%
% Written by EK @ NYU

%% Define all subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'}; %

% Define additional free params / options
addOffsetParam  = false; % Add an offset parameter for each sensor
refitGainParam  = false; % Refit the gain parameter instead of cross-validate

%% Fit model with fMRI pRF parameters
opt = getOpts('saveFig', false,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs', false, ...
    'addOffsetParam', addOffsetParam, 'refitGainParam', refitGainParam); % see getOpts function for more options

for s = length(subjects)
    subjID = subjects{s};
    mprf_main(subjID, opt);
end

%% Fit model with altered pRF size
opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs','size', ...
    'addOffsetParam', addOffsetParam, 'refitGainParam', refitGainParam); % see getOpts function for more options
for s = 1:length(subjects)
    subjID = subjects{s};
    mprf_main(subjID, opt);
    close all
end

%% Fit model with rotated pRF position
opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs','position', ...
    'addOffsetParam', addOffsetParam, 'refitGainParam', refitGainParam); % see getOpts function for more options
for s = 1:length(subjects)
    subjID = subjects{s};
    mprf_main(subjID, opt);
    close all
end


%% Fit model with scrambled pRF parameters across cortex.
%       opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs','scramble', ...
%           'addOffsetParam', addOffsetParam, 'refitGainParam', refitGainParam); % see getOpts function for more options
%       for s = 1:length(subjects)
%           subjID = subjects{s};
%           mprf_main(subjID, opt);
%           close all
%       end


return
