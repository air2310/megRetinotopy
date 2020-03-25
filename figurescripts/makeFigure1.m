function makeFigure1(subjID, dirPth, opt, plotAverage)
% Wrapper function to generate figures 4 A, B and C for the paper.

% INPUTS:
%   opt         : options for modelfit etc (see getOpts)
%   dirPth      : paths with data files for this subject    

% Check if the folder to save the final figures exists, else create one

saveSubDir = ['figure1_' opt.regressionType];

if ~exist(fullfile(dirPth.finalFig.savePth,saveSubDir),'dir')
    mkdir(fullfile(dirPth.finalFig.savePth,saveSubDir));
end

if plotAverage
    
    subjIDs = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058', ...
        'wlsubj068', 'wlsubj070', 'wlsubj081', 'wlsubj106', ...
        'wlsubj109', 'wlsubj111'};
    
%     makeFigure1C(subjIDs,saveSubDir); % Average MEG headplot for all subjects
    mprf_CompareGroupAverageDataVsPrediction(subjIDs, dirPth, opt, saveSubDir)
    
    makeFigure1Supplement(subjIDs,saveSubDir)
    
else
    
    % Plot variance explained by model prediction on a mesh for individual
    % subject
    makeFigure1A(dirPth,opt,saveSubDir);
    
    % Plot 10 MEG time series with pRF prediction for individual subject
    makeFigure1B(dirPth,opt,saveSubDir);
end

end
