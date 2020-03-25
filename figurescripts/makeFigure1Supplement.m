function makeFigure1Supplement(subjIDs,saveSubDir)
% Function to create figure 1A (MEG head plot showing the variance
% explained values for individual subjects).

% define initial parameters
%numChans = length(opt.meg.dataChan);
numChans = 157;
numSub   = length(subjIDs);
varExplAllSubj = nan(numSub,numChans);

maxClim = true;
interpMethod =  'nearest';

% Plot all individuals
fH1 = figure; set(gcf,'Position',[1000, 651, 1500, 687]);

% loop over all subjects
for idxSubj = 1:numSub
    
    % Load paths with data files for this subject
    dirPth = loadPaths(subjIDs{idxSubj});
    opt = getOpts('saveFig',1,'verbose',1, 'fullSizeMesh', 1);
    
    
    varExpFile = dir(fullfile(dirPth.model.saveDataPth,opt.subfolder,'pred_resp','meanVarExpl.mat'));
    if ~isempty(varExpFile)
        load(fullfile(varExpFile.folder,varExpFile.name),'meanVarExpl');       
        varExplAllSubj(idxSubj,:) = meanVarExpl;      
    end
    
    if maxClim
        cmax = max(varExplAllSubj(idxSubj,:));
    else
        cmax = 0.4;
    end
    
    % plot mesh of subject
    subplot(2,5,idxSubj);
    ttl = sprintf('S%d', idxSubj);
    megPlotMap(varExplAllSubj(idxSubj,:),[0 cmax], fH1, 'parula', ttl,[],[],'interpmethod', interpMethod); hold on;
    c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
    pos = c.Position; set(c, 'Position', [pos(1)+0.04 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
    
end

% calculate average variance explained across subjects
meanVarExplAllSubj    = nanmean(varExplAllSubj,1);
[maxMeanVarExplAllSubj, maxSen] = nanmax(meanVarExplAllSubj);

if maxClim
    cmax = max(meanVarExplAllSubj);
else
    cmax = 0.3;
end

% Plot average
fH2 = figure; set(gcf,'Position',[100 100 1920/2 1920/2]);
ttl = strcat('Mean variance explained (Group Average N=', sprintf('%d)',numSub));
megPlotMap(meanVarExplAllSubj,[0 cmax],fH2, 'parula', ttl,[],[],'interpmethod', interpMethod);
c = colorbar; c.TickDirection = 'out'; c.Box = 'off';


if opt.saveFig
    saveDir = fullfile(mprf_rootPath,'data','Retinotopy','Quality_check','average','finalfig',saveSubDir, 'figure1C');
    
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    
    if cmax
        figNameEnd = '_max';
    else
        figNameEnd = '';
    end
        
    figure(fH1);
    figurewrite(fullfile(saveDir, sprintf('fig1Supplement_Individual_subjecs_VE_%s%s', opt.fNamePostFix,figNameEnd)),[],[1 300],'.',1);
    figurewrite(fullfile(saveDir, sprintf('fig1Supplement_Individual_subjecs_VE_%s%s', opt.fNamePostFix,figNameEnd)),[],0,'.',1);
    fprintf('\n saving figure 1D in %s',saveDir);
    
    figure(fH2);
    figurewrite(fullfile(saveDir, sprintf('fig1Supplement_SensorwiseAverage_%d_subjs_maxVar_%2.0f_sen_%d_%s%s',numSub,100*maxMeanVarExplAllSubj,maxSen, opt.fNamePostFix,figNameEnd)),[],[1 300],'.',1);
    figurewrite(fullfile(saveDir, sprintf('fig1Supplement_SensorwiseAverage_%d_subjs_maxVar_%2.0f_sen_%d_%s%s',numSub,100*maxMeanVarExplAllSubj,maxSen, opt.fNamePostFix,figNameEnd)),[],0,'.',1);
    fprintf('\n saving figure 1C in %s',saveDir);
    
    fprintf('\n');
end
end