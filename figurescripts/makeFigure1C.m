function makeFigure1C(subjIDs,opt)
% Function to create figure 1A (MEG head plot showing the variance
% explained values for individual subjects).

% define initial parameters
numChans = length(opt.dataChan);
numSubjs = length(subjIDs);
varExplAllSubj = nan(numSubjs,numChans);

% loop over all subjects
for idxSubj = 1:numSubjs
    % calculate average variance explained across subjects
    
    % Load paths with data files for this subject
    dirPth = loadPaths(subjIDs{idxSubj});
    
    varExpFile = dir(fullfile(dirPth.model.saveDataPth,'original','pred_resp','meanVarExpl.mat'));
    if ~isempty(varExpFile)
        load(fullfile(varExpFile.folder,varExpFile.name),'meanVarExpl');
        
        varExplAllSubj(idxSubj,:) = meanVarExpl;
        
    end
     
end

meanVarExplAllSubj = nanmean(varExplAllSubj,1);

close all;

fH1 = figure; set(gcf,'Position',[100 100 1920/2 1920/2]);
ttl = strcat('Mean variance explained (Average: ', sprintf('%s ,',subjIDs{:}),')');
megPlotMap(meanVarExplAllSubj,[0 0.6],fH1, 'parula', ttl, [],[]);
c = colorbar; c.Location='southoutside';

if opt.saveFig
    saveSubDir = 'figure1C';
    saveDir = fullfile(mprf_rootPath,'data','Retinotopy','Quality_check','averageSub','finalFig',saveSubDir);
    
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    
    print(fH1, fullfile(saveDir, sprintf('Mean_VE_average_%d_subjs',numSubjs)), '-dpng');
    fprintf('\n saving figure 1A in %s',saveDir);
    
    
    fprintf('\n');
end
end