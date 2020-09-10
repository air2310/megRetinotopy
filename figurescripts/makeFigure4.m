function makeFigure4(dirPth, opt)
% Wrapper function to generate figures 4 A, B for the paper.

% INPUTS:
%   opt         : options for modelfit etc (see getOpts)
%   dirPth      : paths with data files for this subject    

% Check if the folder to save the final figures exists, else create one

saveSubDir = ['Figure4_' opt.subfolder];

if ~exist(fullfile(dirPth.finalFig.savePth,saveSubDir),'dir')
    mkdir(fullfile(dirPth.finalFig.savePth,saveSubDir));
end

   
    
    % Plot 10 MEG time series with pRF prediction for individual subject
    makeFigure4A(dirPth,opt,saveSubDir);
    
    % Plot variance explained by model prediction on a mesh for individual
    % subject
    makeFigure4B(dirPth,opt,saveSubDir);
    
end
