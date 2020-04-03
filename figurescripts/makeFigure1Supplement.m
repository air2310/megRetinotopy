function makeFigure1Supplement(dirPth, opt)
% Function to create figure 1A (MEG head plot showing the variance
% explained values for individual subjects).

saveSubDir = ['figure1_' opt.subfolder];

if ~exist(fullfile(dirPth.finalFig.savePth,saveSubDir),'dir')
    mkdir(fullfile(dirPth.finalFig.savePth,saveSubDir));
end


% define initial parameters
subjIDs = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058', ...
    'wlsubj068', 'wlsubj070', 'wlsubj081', 'wlsubj106', ...
    'wlsubj109', 'wlsubj111'};

nChans = 157;
nSub   = length(subjIDs);
varExplAllSubj = nan(nSub,nChans);

% Plotting params
maxClim      = false;
interpMethod =  'nearest';

fH1 = figure; clf; set(gcf,'Position',[1000, 651, 1500, 687]);

% loop over all subjects
for idxSubj = 1:nSub
    
    % Load paths with data files for this subject
    dirPth = loadPaths(subjIDs{idxSubj});
    
    load(fullfile(dirPth.model.saveDataPth,opt.subfolder,'pred_resp','meanVarExpl.mat'), 'meanVarExpl');
    
    varExplAllSubj(idxSubj,:) = meanVarExpl;
end

if maxClim
    cmax = max(varExplAllSubj(:));
else
    cmax = 0.5;
end

for idxSubj = 1:nSub
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
    cmax = max(meanVarExplAllSubj(:));
else
    cmax = 0.3;
end

% Plot average
fH2 = figure; set(gcf,'Position',[100 100 1920/2 1920/2]);
ttl = sprintf('Sensorwise average variance explained (N=%d)',nSub);
megPlotMap(meanVarExplAllSubj,[0 cmax],fH2, 'parula', ttl,[],[],'interpmethod', interpMethod);
c = colorbar; c.TickDirection = 'out'; c.Box = 'off';


if opt.saveFig
    saveDir = fullfile(mprf_rootPath,'data','Retinotopy','Quality_check','average','finalfig',saveSubDir, 'figure1C');
    
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    
    if maxClim
        figNameEnd = '_max';
    else
        figNameEnd = '';
    end
    
    figure(fH1);
    figurewrite(fullfile(saveDir, sprintf('fig1Supplement_Individual_subjecs_VE_%s%s', opt.fNamePostFix,figNameEnd)),[],[1 300],'.',1);
    figurewrite(fullfile(saveDir, sprintf('fig1Supplement_Individual_subjecs_VE_%s%s', opt.fNamePostFix,figNameEnd)),[],0,'.',1);
    fprintf('\n saving supplement to figure 1 in %s',saveDir);
    
    figure(fH2);
    figurewrite(fullfile(saveDir, sprintf('fig1Supplement_SensorwiseAverage_N%d_maxVar_%2.0f_sen_%d_%s%s',nSub,100*maxMeanVarExplAllSubj,maxSen, opt.fNamePostFix,figNameEnd)),[],[1 300],'.',1);
    figurewrite(fullfile(saveDir, sprintf('fig1Supplement_SensorwiseAverage_N%d_maxVar_%2.0f_sen_%d_%s%s',nSub,100*maxMeanVarExplAllSubj,maxSen, opt.fNamePostFix,figNameEnd)),[],0,'.',1);
    fprintf('\n saving sensorwise average figure  in %s',saveDir);
    
    fprintf('\n');
end

return