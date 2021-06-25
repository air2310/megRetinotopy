function makeFigure4B(dirPth,opt,saveSubDir)
% Function to create figure 4B (MEG head plot showing the variance
% explained values for individual subjects).

varexpl = dir(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp','meanVarExpl.mat'));

saveDir = fullfile(dirPth.finalFig.savePth,saveSubDir,'Figure4B');
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
clims = [0 0.5]; % or try clims = [0 max(meanVarExpl)];
interpMethod = 'v4'; % choose 'v4' or 'nearest'
interplim = 'convex';

% Plot it!
fH1 = figure; clf;
megPlotMap(meanVarExpl,clims,fH1, 'parula', ...
    sprintf('Mean variance explained %s',dirPth.subjID), [],[], 'interpmethod', interpMethod, 'mask',interplim);
c = colorbar;
c.Location = 'eastoutside';
c.Box = 'off';
c.TickDirection = 'out';
c.TickLength = [0.010 0.010];
c.FontSize = 12;

if opt.saveFig
    fprintf('\n(%s): Saving Figure 4B in %s\n',mfilename, saveDir);
    figurewrite(fullfile(saveDir, sprintf('Figure4B_Mean_variance_explained_%s_%s_%s', opt.fNamePostFix, interpMethod, dirPth.subjID)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('Figure4B_Mean_variance_explained_%s_%s_%s', opt.fNamePostFix, interpMethod, dirPth.subjID)),[], [1 300],'.',1);
end

end