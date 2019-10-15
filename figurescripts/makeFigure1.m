function makeFigure1(dirPth, opt, plotAverage)
% Warpper function to generate figures 1A and 1B for the paper.

% Check if the folder to save the final figures exists, else create one
saveSubDir = 'figure1';

if ~exist(dirPth.finalFig.savePth,'dir')
    mkdir(dirPth.finalFig.savePth);
end

if ~exist(fullfile(dirPth.finalFig.savePth,saveSubDir),'dir')
    mkdir(fullfile(dirPth.finalFig.savePth,saveSubDir));
end

%%

if plotAverage
    subjIDs = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058', ...
               'wlsubj068', 'wlsubj070', 'wlsubj081', 'wlsubj106', ...
               'wlsubj109', 'wlsubj111'};

    makeFigure1C(subjIDs,opt); % Average MEG headplot for all subjects 
else
    % Plot 10 MEG time series with pRF prediction for individual subject
    makeFigure1A(dirPth,opt);
    
    % Plot variance explained by model prediction on a mesh for individual
    % subject
    makeFigure1B(dirPth,opt);
end

end
