% mprf_runAllAnalysis.m
%
% Master script to run 3 types of analyses in the MEG retinotopy project.
% Default script runs the modelfits with 2 free parameters (beta (or gain)
% and reference phase. Model uses split-half cross-validation to get the 
% the final reference phase and betas within a subject.
%
% Previous versions of the model could add an additional offset, 
% ('addOffsetParam') or refits gain parameters (so only cross-validate
% reference phase ('refitGainParam'). But this model fit can be subject to
% overfitting and is therefore not recommended.

% Define all subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'}; %

% Define additional free params / options
addOffsetParam  = false;
refitGainParam  = false;


% With fMRI pRF parameters
opt = getOpts('saveFig', true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs', false, ...
    'addOffsetParam', addOffsetParam, 'refitGainParam', refitGainParam); % see getOpts function for more options

for s = 1:length(subjects)
    subjID = subjects{s};
    mprf_main(subjID, opt);
end

% With altered pRF size
opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs','size', ...
    'addOffsetParam', addOffsetParam, 'refitGainParam', refitGainParam); % see getOpts function for more options
for s = 1:length(subjects)
    subjID = subjects{s};
    mprf_main(subjID, opt);
    close all
end

% With rotated pRF position
opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs','position', ...
    'addOffsetParam', addOffsetParam, 'refitGainParam', refitGainParam); % see getOpts function for more options
for s = 1:length(subjects)
    subjID = subjects{s};
    mprf_main(subjID, opt);
    close all
end


%       % Scramble
%       opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs','scramble', ...
%           'addOffsetParam', addOffsetParam, 'refitGainParam', refitGainParam); % see getOpts function for more options
%       for s = 1:length(subjects)
%           subjID = subjects{s};
%           mprf_main(subjID, opt);
%           close all
%       end



%% Visualize all types of analyses

% Make individual subject figures
opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs', false,...
    'addOffsetParam', addOffsetParam, 'refitGainParam', refitGainParam); % see getOpts function for more options

for s = 1:length(subjects)
    subjID = subjects{s};
    close all; makeAllFigures(subjID, 1:3, 'top10', false, 'meanVE', opt);
end

% Make average figures
close all; makeAllFigures('wlsubj004', 1:3, 'top10', true, 'meanVE', opt);


return
