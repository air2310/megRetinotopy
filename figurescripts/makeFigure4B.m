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
if opt.mri.useHCPAveMaps
    clims = [0 max(meanVarExpl,[],'omitnan')];
else
    clims = [0 0.5];
end
% clims = [0 max(meanVarExpl)];
interpMethod = 'v4'; % choose 'v4' or 'nearest'

% Plot it!
fH1 = figure; clf;
megPlotMap(meanVarExpl,clims,fH1, 'parula', ...
    sprintf('Mean variance explained %s',dirPth.subjID), [],[], 'interpmethod', interpMethod);
c = colorbar;
c.Location = 'eastoutside';
c.Box = 'off';
c.TickDirection = 'out';
c.TickLength = [0.010 0.010];
c.FontSize = 12;

if opt.saveFig
    fprintf('\n(%s): Saving figure 4B in %s\n',mfilename, saveDir);
    figurewrite(fullfile(saveDir, sprintf('Figure4B_Mean_variance_explained_%s_%s', opt.fNamePostFix, interpMethod)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('Figure4B_Mean_variance_explained_%s_%s', opt.fNamePostFix, interpMethod)),[], [1 300],'.',1);
end

end