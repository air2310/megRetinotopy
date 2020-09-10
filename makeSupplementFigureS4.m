function makeSupplementFigureS4
% Function to visualize meshes for group average fit when rotating pRF position

% Written by Eline Kupers, NYU 2020

opt = getOpts('verbose', verbose, 'saveFig', saveFig, ...
    'perturbOrigPRFs', 'position', 'headmodel', headmodel, ...
    'fullSizeMesh',fullSizeMesh, 'addOffsetParam', addOffsetParam, ...
    'refitGainParam', refitGainParam);

%% Make save figure dir
saveSubDir = sprintf('SupplFigureS4');
saveDir = fullfile(dirPth.finalFig.savePthAverage,saveSubDir);
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

%% Handle data

% Define all subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

% Load data
[dataDir, tmp] = fileparts(dirPth.finalFig.savePthAverage);
loadDir = fullfile(dataDir, 'GroupAvePrediction', opt.subfolder);

% Check number of iterations
if ~opt.vary.perturbOrigPRFs
    nrVariations = 1;
elseif strcmp(opt.vary.perturbOrigPRFs, 'position')
    nrVariations = length(opt.vary.position);
elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
    nrVariations = length(opt.vary.size);
end

nSensors = 157;
nBoot = 10000;
groupVarExplAll = NaN(nrVariations, nSensors, nBoot);
groupVarExplMean = NaN(nrVariations, nSensors);


for ii = 1:nrVariations
    load(fullfile(loadDir, sprintf('groupVarExplBoot10000_%d',ii)),'groupVarExpl', 'groupAveData', 'groupAvePredictionScaled', 'opt', 'nBoot');
    load(fullfile(loadDir, sprintf('groupVarExplNoBoot_%d',ii)),'groupVE_noBootstrp', 'groupAvePredScaled_noBootstrp', 'allData', 'allPredictions');
    
    if ~opt.vary.perturbOrigPRFs
        groupVarExplMeanBoot = nanmean(groupVarExpl,2);
        ve_cat = cat(1,groupVarExplMeanBoot',groupVE_noBootstrp);
    else
        groupVarExplAll(ii,:,:) = groupVarExpl;
        groupVarExplMean(ii,:) = nanmean(groupVarExpl,2)';
    end
end


% Plot var explained map
clims = [0 0.5];
fH1 = figure(1); clf; drawnow; set(fH1, 'Position', [1000,420,1269,918]);
fH2 = figure(2); clf; drawnow; set(fH2, 'Position', [1000,420,1269,918]);

for ii = 1:nrVariations
    
    figure(fH1); subplot(3,ceil(nrVariations/3),ii)
    megPlotMap(groupVarExplMean(ii,:),clims,fH1, 'parula', ...
        subplotLabels{ii}, [],[], 'interpMethod', 'v4');
    c = colorbar;
    c.Location = 'eastoutside';
    c.Box = 'off';
    c.TickDirection = 'out';
    c.TickLength = [0.010 0.010];
    c.FontSize = 12;
    pos = c.Position; set(c, 'Position', [pos(1)+0.03 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
    
%     figure(fH2); subplot(3,ceil(nrVariations/3),ii)
%     megPlotMap(groupVarExplMean(ii,:),clims,fH2, 'parula', ...
%         subplotLabels{ii}, [],[], 'interpmethod', 'nearest');
%     c = colorbar;
%     c.Location = 'eastoutside';
%     c.Box = 'off';
%     c.TickDirection = 'out';
%     c.TickLength = [0.010 0.010];
%     c.FontSize = 12;
%     pos = c.Position; set(c, 'Position', [pos(1)+0.03 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
%     
end

% figure(fH2)
% figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_vary%s_nearest',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs)),[],0,'.',1);
% figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_vary%s_nearest',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs)),[],[1 300],'.',1);

figure(fH1)
figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_vary%s',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs)),[],0,'.',1);
figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_vary%s',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs)),[],[1 300],'.',1);


return