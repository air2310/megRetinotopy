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
makeFigure1A(dirPth);

%%
makeFigure1B(dirPth,opt);


close all;
end