function makeFigure1A(dirPth,opt)
% Function to create figure 1A (MEG head plot showing the variance
% explained values for individual subjects).
varexpl = dir(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp','meanVarExpl.mat'));


saveSubDir = 'figure1A';
saveDir = fullfile(dirPth.finalFig.savePth,'figure1',saveSubDir);
if ~exist(saveDir,'dir')
    mkdir(fullfile(dirPth.finalFig.savePth,'figure1',saveSubDir));
end

% check if the modelPredictions are saved in the folder.
close all;
if isempty(varexpl.folder)
    error('Can''t find mean variance explained mat-file')
end

% Load the variance explained file
load(fullfile(varexpl.folder,varexpl.name),'meanVarExpl');

% Set colormap limits
clims = [0 max(meanVarExpl)];

% Plot it!
fH1 = figure; clf;
megPlotMap(meanVarExpl,clims,fH1, 'parula', ...
    'Mean variance explained', [],[], 'interpmethod', 'nearest');
c = colorbar;
c.Location = 'southoutside';
c.Box = 'off';
c.TickDirection = 'out';
c.TickLength = [0.010 0.010];
c.FontSize = 12;


if opt.saveFig
    print(fH1, fullfile(saveDir, sprintf('Mean_variance_explained_%s', opt.fNamePostFix)), '-dpng');
    saveas(fH1, fullfile(saveDir, sprintf('Mean_variance_explained_%s', opt.fNamePostFix')), 'epsc');
    fprintf('\n(%s): Saving figure 1A in %s\n',mfilename, saveDir);
end

end