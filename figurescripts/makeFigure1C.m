function makeFigure1C(subjIDs)
% Function to create figure 1A (MEG head plot showing the variance
% explained values for individual subjects).

% define initial parameters
numChans = 157;
numSub = length(subjIDs);
varExplAllSubj = nan(numSub,numChans);

% loop over all subjects
for idxSubj = 1:numSub
    % calculate average variance explained across subjects
    
    % Load paths with data files for this subject
    dirPth = loadPaths(subjIDs{idxSubj});
    opt = getOpts('saveFig',1,'verbose',1, 'fullSizeMesh', 1);
    
    
    varExpFile = dir(fullfile(dirPth.model.saveDataPth,opt.subfolder,'pred_resp','meanVarExpl.mat'));
    if ~isempty(varExpFile)
        load(fullfile(varExpFile.folder,varExpFile.name),'meanVarExpl');
        
        varExplAllSubj(idxSubj,:) = meanVarExpl;
        
    end
     
end

meanVarExplAllSubj    = nanmean(varExplAllSubj,1);
[maxMeanVarExplAllSubj, maxSen] = nanmax(meanVarExplAllSubj);

close all;


% Individual subjects MEG head plots
ttl = strcat('Variance explained ( subjects: N=', sprintf('%d )',numSub));
fH1 = figure; set(gcf,'Position',[66,1,1855,1001],'Name',ttl);
for i = 1:numSub
    subplot(2,5,i);
    megPlotMap(varExplAllSubj(i,:),[0 maxMeanVarExplAllSubj],fH1, 'parula',[],[],[],'interpmethod', 'nearest');
    c = colorbar; c.Location='southoutside';
end


% Average across subjects
fH2 = figure; set(gcf,'Position',[100 100 1920/2 1920/2]);
ttl = strcat('Mean variance explained (Average across subjects: N=', sprintf('%d )',numSub));
megPlotMap(meanVarExplAllSubj,[0 maxMeanVarExplAllSubj],fH2, 'parula', ttl,[],[],'interpmethod', 'nearest');
c = colorbar; c.Location='southoutside';


if opt.saveFig
    saveSubDir = 'figure1C';
    saveDir = fullfile(mprf_rootPath,'data','Retinotopy','Quality_check','averageSub','finalFig',saveSubDir);
    
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end

    print(fH1, fullfile(saveDir, sprintf('VE_individual_%d_subjs_maxVar_%d_sen_%d',numSub,round(maxMeanVarExplAllSubj),maxSen)), '-dpng');
    
    print(fH2, fullfile(saveDir, sprintf('Mean_VE_average_%d_subjs_maxVar_%d_sen_%d',numSub,round(maxMeanVarExplAllSubj),maxSen)), '-dpng');
    print(fH2, fullfile(saveDir, sprintf('Mean_VE_average_%d_subjs_maxVar_%d_sen_%d',numSub,round(maxMeanVarExplAllSubj),maxSen)), '-depsc');
    fprintf('\n saving figure 1C in %s',saveDir);
    
    
    fprintf('\n');
end
end