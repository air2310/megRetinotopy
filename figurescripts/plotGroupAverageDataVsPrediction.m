function plotGroupAverageDataVsPrediction(dirPth, opt, whichFigure, sensorsToAverage)
% Function to visualize group average fit (an alternative for the
% sensor-wise average across individual subjects)
%
% This function depends on running the actual analysis, which can be done
% for the different analysis types setting opt.
% Example 1, using original pRFs:
%   opt = getOpts('perturbOrigPRFs',false);
%   mprf_CompareGroupAverageDataVsPrediction(opt);
%
% Example 2, with pRF position change:
%   opt = getOpts('perturbOrigPRFs','position');
%   mprf_CompareGroupAverageDataVsPrediction(opt);
%
% To run this plotting function use the wrapper 'makeAllFigures', e.g.:
%   makeAllFigures('wlsubj004', 1, 'top10', true, 'meanVE')
%
%
% Written by Eline Kupers, NYU 2020

%% Check inputs
if ~exist('sensorsToAverage', 'var') || isempty(sensorsToAverage)
    sensorsToAverage = 'top10';
end

if ~exist('plotSupplementFig', 'var') || isempty(plotSupplementFig)
    plotSupplementFig = false;
end


%% Make save figure dir
saveSubDir = sprintf('Figure%d_%s',whichFigure, opt.subfolder);
saveDir = fullfile(dirPth.finalFig.savePthAverage,saveSubDir, 'GroupAvePrediction');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

%% Handle data

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
interpMthd = 'v4';

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

if ~opt.vary.perturbOrigPRFs
    fnameAdd = {'withBoot', 'noBoot'};
    for ii = [1,2]
        groupVarExplMean =  ve_cat(ii,:);
        
        % Plot var explained map
        clims = [0 max(groupVarExplMean)];
        fH1 = figure; clf;
        megPlotMap(groupVarExplMean,clims,fH1, 'parula', ...
            sprintf('Average Group Fit Variance Explained %s', fnameAdd{ii}), [],[], 'interpMethod', interpMthd);
        c = colorbar;
        c.Location = 'eastoutside';
        c.Box = 'off';
        c.TickDirection = 'out';
        c.TickLength = [0.010 0.010];
        c.FontSize = 12;
        
        figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_%s_%s',opt.fNamePostFix, opt.addOffsetParam, fnameAdd{ii}, interpMthd)),[],0,'.',1);
        figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_%s_%s',opt.fNamePostFix, opt.addOffsetParam, fnameAdd{ii}, interpMthd)),[],[1 300],'.',1);
        
    end
    
else
    
    if strcmp(opt.vary.perturbOrigPRFs,'position')
        range = opt.vary.position;
        origIdx = find(range==0);
        for l = 1:length(range)
            subplotLabels{l} = sprintf('%d deg', rad2deg(range(l)));
        end
        xLabels = 'Rotation angle (deg)';
        xScale = 'linear';
    elseif strcmp(opt.vary.perturbOrigPRFs,'size')
        range = opt.vary.size;
        origIdx = find(range==1);
        for l = 1:length(range)
            subplotLabels{l} = sprintf('%1.1fx', opt.vary.size(l));
        end
        xLabels = 'Size Scale factor';
        xScale = 'log';
    end
    
    
    % Select sensors to average
    sensorLoc = selectSensorsToAverage(opt, dirPth, saveDir, groupVarExplMean, sensorsToAverage);
    
    % Select data to plot
    if strcmp(sensorsToAverage, 'top10Positive')
        curSubjLocs = sensorLoc;
        for ii = 1:size(curSubjLocs,1)
            curSensorLoc = curSubjLocs(ii,:);
            dataToPlot(ii,1:sum(~isnan(curSensorLoc))) = squeeze(groupVarExplMean(curSensorLoc(~isnan(curSensorLoc)),ii));
        end
    else
        dataToPlot = squeeze(groupVarExplMean(:,sensorLoc))';
    end
    
    % Get confidence intervals based on boots across subjects
    pct1 = 100 * (0.32/2);
    pct2 = 100 - pct1;
    lo = 100.*nanmean(prctile(groupVarExplAll(:,sensorLoc,:), pct1, 3),2); % 16th percentile
    hi = 100.*nanmean(prctile(groupVarExplAll(:,sensorLoc,:), pct2, 3),2); % 84th percentile
    
    % Get p value from bootstrapped data
    mnTop10Boots = squeeze(nanmean(groupVarExplAll(:,sensorLoc,:),2));
    [~,peakIdx] = max(mnTop10Boots,[],1);
    peakBelowOrigPRF = sum(peakIdx<origIdx);
    p = 2*(.5-abs(.5- (peakBelowOrigPRF/size(groupVarExplAll,3))));
    fprintf('P-value of mean peak at original pRF param using %dx bootstraps: %1.3f \n',size(groupVarExplAll,3), p)
    
    
    figure(99); set(gcf, 'Position',[520, 39, 996, 759], 'Color', 'w'); clf;
    meandataToPlot = 100.*nanmean(dataToPlot, 1);
    plot([range(origIdx) range(origIdx)],[0, max(hi')], 'k');  hold on;
    plot([range(1), range(end)], [0 0], 'k');
    err = patch([range, fliplr(range)], [lo', fliplr(hi')], [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'LineStyle',':');
    plot(range,meandataToPlot,'Color', 'r', 'Linewidth',2); hold on;
    
    set(gca,'TickDir', 'out', 'YLim', [0, max(hi)], 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20);
    xlabel(xLabels);
    set(gca,'XScale', xScale, 'XLim', [range(1),range(end)]); axis square;
    title('Group average fit');
    ylabel('Variance explained (%)'); box off;
    
    figurewrite(fullfile(saveDir, sensorsToAverage, sprintf('GroupAverageFit_VarExplLine_%s_Offset%d_vary%s_%s',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs,sensorsToAverage)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sensorsToAverage, sprintf('GroupAverageFit_VarExplLine_%s_Offset%d_vary%s_%s',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs,sensorsToAverage)),[],[1 300],'.',1);
    
end

    
return