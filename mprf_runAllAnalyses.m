% mprf_runAllAnalysis.m
% 
% Master script to run 4 types of analyses in the MEG retinotopy project.
% The script runs the modelfits with either 2 free parameters (reference
% phase and gain, 'NoOffset') or 3 free parameters (reference phase, gain
% and offset, 'WithOffset')

% Define all subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'}; %

addOffsetParam  = [false, true];
refitGainParam  = [false, true];

for rt = 1:length(addOffsetParam)
      offset = addOffsetParam(rt);
    
    for rb = 1:length(refitGainParam)
        refitGain = refitGainParam(rb);
        
        
        % With fMRI pRF parameters
        opt = getOpts('saveFig', true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs', false, ...
            'addOffsetParam', offset, 'refitGainParam', refitGain); % see getOpts function for more options          
        
        for s = 1:length(subjects)
            subjID = subjects{s};
            mprf_main(subjID, opt);
        end
        
        % With altered pRF size
        opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs','size', ...
            'addOffsetParam', offset, 'refitGainParam', refitGain); % see getOpts function for more options
        for s = 1:length(subjects)
            subjID = subjects{s};
            mprf_main(subjID, opt);
            close all
        end
        
        % With rorated pRF position
        opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs','position', ...
            'addOffsetParam', offset, 'refitGainParam', refitGain); % see getOpts function for more options
        for s = 1:length(subjects)
            subjID = subjects{s};
            mprf_main(subjID, opt);
            close all
        end
        
        
%       % Scramble
%       opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs','scramble', ...
%           'addOffsetParam', addOffsetParam, 'refitGainParam', refitGain); % see getOpts function for more options
%       for s = 1:length(subjects)
%           subjID = subjects{s};
%           mprf_main(subjID, opt);
%           close all
%       end
        
    end
end


%% Visualize all types of analyses
for rt = 1:length(addOffsetParam)
    offset = addOffsetParam(rt);
    
    for rb = 1:length(refitGainParam)
        refitGain = refitGainParam(rb);
        
        
        % Make individual subject figures
        opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs', false,...
            'addOffsetParam', offset, 'refitGainParam', refitGain); % see getOpts function for more options
        
        for s = 1:length(subjects)
            subjID = subjects{s};
            close all; makeAllFigures(subjID, 1:3, 'top10', false, 'meanVE', opt);
        end
        
        % Make average figures
        close all; makeAllFigures('wlsubj004', 1:4, 'top10', true, 'meanVE', opt);
    end
end

return
