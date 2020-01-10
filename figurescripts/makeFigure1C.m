function makeFigure1C(subjIDs)
% Function to create figure 1A (MEG head plot showing the variance
% explained values for individual subjects).

% define initial parameters
%numChans = length(opt.meg.dataChan);
numChans = 157;
numSub   = length(subjIDs);
varExplAllSubj = nan(numSub,numChans);

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
    
    % plot mesh of subject
    subplot(2,5,idxSubj);
    ttl = sprintf('S%d', idxSubj);
    megPlotMap(varExplAllSubj(idxSubj,:),[0 0.45], fH1, 'parula', ttl,[],[],'interpmethod', []); hold on;
    c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
    pos = c.Position; set(c, 'Position', [pos(1)+0.04 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
    
end

% calculate average variance explained across subjects
meanVarExplAllSubj    = nanmean(varExplAllSubj,1);
[maxMeanVarExplAllSubj, maxSen] = nanmax(meanVarExplAllSubj);

% Plot average
fH2 = figure; set(gcf,'Position',[100 100 1920/2 1920/2]);
ttl = strcat('Mean variance explained (Group Average N=', sprintf('%d)',numSub));
megPlotMap(meanVarExplAllSubj,[0 0.30],fH2, 'parula', ttl,[],[],'interpmethod', []);
c = colorbar; c.TickDirection = 'out'; c.Box = 'off';


if opt.saveFig
    saveSubDir = 'figure1C';
    saveDir = fullfile(mprf_rootPath,'data','Retinotopy','Quality_check','average','finalfig',saveSubDir);
    
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    
    figure(fH1);
    figurewrite(fullfile(saveDir, sprintf('fig1D_Individual_subjecs_VE_%s', opt.fNamePostFix)),[],[1 300],'.',1);
    figurewrite(fullfile(saveDir, sprintf('fig1D_Individual_subjecs_VE_%s', opt.fNamePostFix)),[],0,'.',1);
    fprintf('\n saving figure 1D in %s',saveDir);
    
    figure(fH2);
    figurewrite(fullfile(saveDir, sprintf('fig1C_Mean_VE_average_%d_subjs_maxVar_%2.0f_sen_%d_%s',numSub,100*maxMeanVarExplAllSubj,maxSen, opt.fNamePostFix)),[],[1 300],'.',1);
    figurewrite(fullfile(saveDir, sprintf('fig1C_Mean_VE_average_%d_subjs_maxVar_%2.0f_sen_%d_%s',numSub,100*maxMeanVarExplAllSubj,maxSen, opt.fNamePostFix)),[],0,'.',1);
    fprintf('\n saving figure 1C in %s',saveDir);
    
    fprintf('\n');
end
end