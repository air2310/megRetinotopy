function makeFigure1A(dirPth,opt,saveSubDir)
% Function to create figure 1A (MEG head plot showing the variance
% explained values for individual subjects).

varexpl = dir(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp','meanVarExpl.mat'));

saveDir = fullfile(dirPth.finalFig.savePth,saveSubDir,'figure1A');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% check if the modelPredictions are saved in the folder.
close all;
if isempty(varexpl.folder)
    error('Can''t find mean variance explained mat-file')
end

% Load the variance explained file
load(fullfile(varexpl.folder,varexpl.name),'meanVarExpl');

% Set colormap limits
clims = [0 0.45];
clims = [0 max(meanVarExpl)];
interpMethod = 'v4'; % choose 'v4' or 'nearest'

% Plot it!
fH1 = figure; clf;
megPlotMap(meanVarExpl,clims,fH1, 'parula', ...
    'Mean variance explained', [],[], 'interpmethod', interpMethod);
c = colorbar;
c.Location = 'eastoutside';
c.Box = 'off';
c.TickDirection = 'out';
c.TickLength = [0.010 0.010];
c.FontSize = 12;

if opt.saveFig
    fprintf('\n(%s): Saving figure 1A in %s\n',mfilename, saveDir);
    figurewrite(fullfile(saveDir, sprintf('Mean_variance_explained_%s_%s', opt.fNamePostFix, interpMethod)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('Mean_variance_explained_%s_%s', opt.fNamePostFix, interpMethod)),[], [1 300],'.',1);
end

end