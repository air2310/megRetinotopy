function plotpRFPerturbationMeshesGroupAverageFit(dirPth, opt)
% Function to make third row in Supplemental Figure S$ or S6 of manuscript,
% showing the effect of perturbing the original pRFs on a head plot, for
% the group average modelfit


%% Load data
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


if strcmp(opt.vary.perturbOrigPRFs,'position')
    range = opt.vary.position;
    origIdx = find(range==0);
    for l = 1:length(range)
        subplotLabels{l} = sprintf('%d deg', rad2deg(range(l)));
    end
    xLabels = 'Rotation angle (deg)';
    xScale = 'linear';
    figNum = 4;
elseif strcmp(opt.vary.perturbOrigPRFs,'size')
    range = opt.vary.size;
    origIdx = find(range==1);
    for l = 1:length(range)
        subplotLabels{l} = sprintf('%1.1fx', opt.vary.size(l));
    end
    xLabels = 'Size Scale factor';
    xScale = 'log';
    figNum = 6;
end


% Plot var explained map
clims = [0 0.5];
interpMthd = 'v4'; % or 'nearest';
fH1 = figure(1); clf; drawnow; set(fH1, 'Position', [1000,420,1269,918]);

for ii = 1:nrVariations
    
    figure(fH1); subplot(3,ceil(nrVariations/3),ii)
    megPlotMap(groupVarExplMean(ii,:),clims,fH1, 'parula', ...
        subplotLabels{ii}, [],[], 'interpMethod', interpMthd);
    c = colorbar;
    c.Location = 'eastoutside';
    c.Box = 'off';
    c.TickDirection = 'out';
    c.TickLength = [0.010 0.010];
    c.FontSize = 12;
    pos = c.Position; set(c, 'Position', [pos(1)+0.03 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
    
end

%% Make save figure dir
if opt.saveFig
    saveSubDir = sprintf('SupplFigureS%d_%s', figNum, opt.subfolder);
    saveDir = fullfile(dirPth.finalFig.savePthAverage,saveSubDir, 'GroupAvePrediction');
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    
    figure(fH2)
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_vary%s_%s',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs, interpMthd)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_vary%s_%s',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs, interpMthd)),[],[1 300],'.',1);
    
end
end
