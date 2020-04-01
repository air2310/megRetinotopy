function makeFigure1(dirPth, opt)
% Wrapper function to generate figures 4 A, B and C for the paper.

% INPUTS:
%   opt         : options for modelfit etc (see getOpts)
%   dirPth      : paths with data files for this subject    

% Check if the folder to save the final figures exists, else create one

saveSubDir = ['figure1_' opt.regressionType];

if ~exist(fullfile(dirPth.finalFig.savePth,saveSubDir),'dir')
    mkdir(fullfile(dirPth.finalFig.savePth,saveSubDir));
end

    
    % Plot variance explained by model prediction on a mesh for individual
    % subject
    makeFigure1A(dirPth,opt,saveSubDir);
    
    % Plot 10 MEG time series with pRF prediction for individual subject
    makeFigure1B(dirPth,opt,saveSubDir);
    
end
