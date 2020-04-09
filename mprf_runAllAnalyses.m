% Define all subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'}; %

% reset random number generator seed for reproducibility
rng('default') % reset to default is sometimes needed before setting specific
% rng(1)

regressionTypes = {'NoOffset', 'WithOffset'};
refitBetas      = [false, true];

for rt = 1:length(regressionTypes)
      regressionType = regressionTypes{rt};
    
    for rb = 1:length(refitBetas)
        recomputeBetaFlag = refitBetas(rb);
        
        
        % With fMRI pRF parameters
        opt = getOpts('saveFig', true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs', false, ...
            'regressionType', regressionType, 'recomputeFinalPredictionBetas', recomputeBetaFlag); % see getOpts function for more options          
        
        for s = 1:length(subjects)
            subjID = subjects{s};
            mprf_main(subjID, opt);
        end
        
        % With altered pRF size
        opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs','size', ...
            'regressionType', regressionType, 'recomputeFinalPredictionBetas', recomputeBetaFlag); % see getOpts function for more options
        for s = 1:length(subjects)
            subjID = subjects{s};
            mprf_main(subjID, opt);
            close all
        end
        
        % With rorated pRF position
        opt = getOpts('saveFig',true,'verbose',false, 'fullSizeMesh', true, 'perturbOrigPRFs','position', ...
            'regressionType', regressionType, 'recomputeFinalPredictionBetas', recomputeBetaFlag); % see getOpts function for more options
        for s = 1:length(subjects)
            subjID = subjects{s};
            mprf_main(subjID, opt);
            close all
        end
        
        
%       % Scramble
%       opt = getOpts('saveFig',1,'verbose',0, 'fullSizeMesh', 1, 'perturbOrigPRFs','scramble', ...
%           'regressionType', regressionType, 'recomputeFinalPredictionBetas', recomputeBetaFlag); % see getOpts function for more options
%       for s = 1:length(subjects)
%           subjID = subjects{s};
%           mprf_main(subjID, opt);
%           close all
%       end
        
    end
end


%% Visualize all types of analyses
for rt = 1:length(regressionTypes)
    regressionType = regressionTypes{rt};
    
    for rb = 1:length(refitBetas)
        recomputeBetaFlag = refitBetas(rb);
        
        
%         Make figures
        opt = getOpts('saveFig',1,'verbose',0, 'fullSizeMesh', 1, 'perturbOrigPRFs', false,...
            'regressionType', regressionType, 'recomputeFinalPredictionBetas', recomputeBetaFlag); % see getOpts function for more options
        
        for s = 1:length(subjects)
            subjID = subjects{s};
            close all; makeAllFigures(subjID, 1:3, 'top10', false, 'meanVE', opt);
        end
        
%         Make average figure
        close all; makeAllFigures('wlsubj004', 1:3, 'top10', true, 'meanVE', opt);
    end
end

return
