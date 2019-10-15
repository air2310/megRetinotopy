function makeFigure1(dirPth,opt)
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
makeFigure1A(dirPth,opt);

if opt.makeAverageFig
    subjIDs = {'wlsubj004','wlsubj039','wlsubj040','wlsubj058','wlsubj068','wlsubj070','wlsubj081','wlsubj106','wlsubj109','wlsubj111'};
    opt_all.dataChan = 1:157;          % MEG Sensor info
    opt_all.saveFig  = true;
    makeFigure1C(subjIDs,opt_all); % Average MEG headplot for all subjects 
end

%%
makeFigure1B(dirPth,opt);


close all;
end