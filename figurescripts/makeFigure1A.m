function makeFigure1A(dirPth,opt)
% Function to create figure 1A (MEG head plot showing the variance
% explained values for individual subjects).

varExpFile = dir(fullfile(dirPth.model.saveDataPth,'original','coherent','pred_resp','meanVarExpl.mat'));

saveSubDir = 'figure1A';
saveDir = fullfile(dirPth.finalFig.savePth,'figure1',saveSubDir);
if ~exist(saveDir,'dir')
    mkdir(fullfile(dirPth.finalFig.savePth,'figure1',saveSubDir));
end

% check if the modelPredictions are saved in the folder.
close all;
if isempty(varExpFile.folder)
    error('Can''t find mean variance explained mat-file')
end

% Load the variance explained file
load(fullfile(varExpFile.folder,varExpFile.name),'meanVarExpl');

% Set colormap limits
clims = [0 0.4];

% Plot it!
fH1 = figure;
megPlotMap(meanVarExpl,clims,fH1, 'parula', ...
    'mean variance explained', [],[], 'interpmethod', 'nearest');
c = colorbar; c.Location='southoutside';


if opt.saveFig
    print(fH1, fullfile(saveDir, sprintf('Mean_variance_explained')), '-dpng');
    saveas(fH1, fullfile(saveDir, sprintf('Mean_variance_explained')), 'epsc');
    fprintf('\n(%s): Saving figure 1A in %s\n',mfilename, saveDir);
end

end