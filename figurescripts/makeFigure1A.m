function makeFigure1A(dirPth)
% Function to create figure 1A (MEG head plot showing the variance
% explained values)

varExpFile = dir(fullfile(dirPth.model.saveDataPth,'original','pred_resp','meanVarExpl.mat'));

saveSubDir = 'figure1A';
saveDir = fullfile(dirPth.finalFig.savePth,'figure1',saveSubDir);
if ~exist(saveSubDir,'dir')
    mkdir(fullfile(dirPth.finalFig.savePth,'figure1',saveSubDir));
end

% check if the modelPredictions are saved in the folder. Else run
% mprf_main.m
close all;
if ~isempty(varExpFile)
    load(fullfile(varExpFile.folder,varExpFile.name),'meanVarExpl');
    
    fH1 = figure; megPlotMap(meanVarExpl,[0 0.6],fH1, 'parula', 'mean variance explained', [],[]);
    c = colorbar; c.Location='southoutside';    
   
    print(fH1, fullfile(saveDir, sprintf('Mean_variance_explained')), '-dpng');
end


end