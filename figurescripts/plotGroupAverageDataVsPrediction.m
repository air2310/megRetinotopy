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

%% Make save figure dir
saveSubDir = sprintf('figure%d_%s',whichFigure, opt.subfolder);
saveDir = fullfile(dirPth.finalFig.savePthAverage,saveSubDir, 'GroupAvePrediction');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

%% Handle data

% Load data
[dataDir, tmp] = fileparts(dirPth.finalFig.savePthAverage);
loadDir = fullfile(dataDir, 'GroupAvePrediction', opt.subfolder);

load(fullfile(loadDir, 'groupVarExplBoot10000'),'groupVarExpl', 'groupAveData', 'groupAvePredictionScaled', 'opt', 'nBoot');
load(fullfile(loadDir, 'groupVarExplNoBoot'),'groupVE_noBootstrp', 'groupAvePredScaled_noBootstrp', 'allData', 'allPredictions');
    

% Check number of iterations
if ~opt.vary.perturbOrigPRFs
    nrVariations = 1;
    groupVarExplMeanBoot = nanmean(groupVarExpl,2);
    ve_cat = cat(1,groupVarExplMeanBoot',groupVE_noBootstrp);
elseif strcmp(opt.vary.perturbOrigPRFs, 'position')
    nrVariations = length(opt.vary.position);
    groupVarExplMean = nanmean(groupVarExpl,3);

elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
    nrVariations = length(opt.vary.size);
    groupVarExplMean = nanmean(groupVarExpl,3);

end


if ~opt.vary.perturbOrigPRFs
    fnameAdd = {'withBoot', 'noBoot'};
    for ii = [1,2]
        groupVarExplMean =  ve_cat(ii,:);
        
        % Plot var explained map
        clims = [0 max(groupVarExplMean)];
        fH1 = figure; clf;
        megPlotMap(groupVarExplMean,clims,fH1, 'parula', ...
            sprintf('Average Group Fit Variance Explained %s', fnameAdd{ii}), [],[], 'interpMethod', 'v4');
        c = colorbar;
        c.Location = 'eastoutside';
        c.Box = 'off';
        c.TickDirection = 'out';
        c.TickLength = [0.010 0.010];
        c.FontSize = 12;

        figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_%s',opt.fNamePostFix, opt.addOffsetParam, fnameAdd{ii})),[],0,'.',1);
        figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_%s',opt.fNamePostFix, opt.addOffsetParam, fnameAdd{ii})),[],[1 300],'.',1);

        fH1 = figure; clf;
        megPlotMap(groupVarExplMean,clims,fH1, 'parula', ...
            sprintf('Average Group Fit Variance Explained %s', fnameAdd{ii}), [],[], 'interpmethod', 'nearest');
        c = colorbar;
        c.Location = 'eastoutside';
        c.Box = 'off';
        c.TickDirection = 'out';
        c.TickLength = [0.010 0.010];
        c.FontSize = 12;

        figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_nearest_%s',opt.fNamePostFix, opt.addOffsetParam, fnameAdd{ii})),[],0,'.',1);
        figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_nearest_%s',opt.fNamePostFix, opt.addOffsetParam, fnameAdd{ii})),[],[1 300],'.',1);
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
    lo = 100.*nanmean(prctile(groupVarExpl(:,sensorLoc,:), pct1, 3),2); % 16th percentile
    hi = 100.*nanmean(prctile(groupVarExpl(:,sensorLoc,:), pct2, 3),2); % 84th percentile
    
    figure(99); clf;
    meandataToPlot = 100.*nanmean(dataToPlot, 1);   
    plot([range(origIdx) range(origIdx)],[min(lo'), max(hi')], 'k');  hold on;
    plot([range(1), range(end)], [0 0], 'k');
    err = patch([range, fliplr(range)], [lo', fliplr(hi')], [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'LineStyle',':');
    plot(range,meandataToPlot,'Color', 'r', 'Linewidth',2); hold on;
    
    set(gca,'TickDir', 'out', 'YLim', [min(lo), max(hi)], 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20);
    xlabel(xLabels);
    set(gca,'XScale', xScale, 'XLim', [range(1),range(end)]); axis square;
    title('Group average fit');
    ylabel('Variance explained (%)'); box off;
    
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplLine_%s_Offset%d_vary%s_%s',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs,sensorsToAverage)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplLine_%s_Offset%d_vary%s_%s',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs,sensorsToAverage)),[],[1 300],'.',1);
    
    
    % Plot var explained map
    clims = [0 0.5];
    fH1 = figure(1); clf; drawnow; set(fH1, 'Position', [1000,420,1269,918]);
    fH2 = figure(2); clf; drawnow; set(fH2, 'Position', [1000,420,1269,918]);
    
    for ii = 1:nrVariations
        
        figure(fH1); subplot(3,ceil(nrVariations/3),ii)
        megPlotMap(groupVarExplMean(:,ii)',clims,fH1, 'parula', ...
            subplotLabels{ii}, [],[], 'interpMethod', 'v4');
        c = colorbar;
        c.Location = 'eastoutside';
        c.Box = 'off';
        c.TickDirection = 'out';
        c.TickLength = [0.010 0.010];
        c.FontSize = 12;
        pos = c.Position; set(c, 'Position', [pos(1)+0.03 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
        
        figure(fH2); subplot(3,ceil(nrVariations/3),ii)
        megPlotMap(groupVarExplMean(:,ii)',clims,fH2, 'parula', ...
            subplotLabels{ii}, [],[], 'interpmethod', 'nearest');
        c = colorbar;
        c.Location = 'eastoutside';
        c.Box = 'off';
        c.TickDirection = 'out';
        c.TickLength = [0.010 0.010];
        c.FontSize = 12;
        pos = c.Position; set(c, 'Position', [pos(1)+0.03 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
        
    end
    
    figure(fH1)
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_vary%s',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_vary%s',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs)),[],[1 300],'.',1);
    
    figure(fH2)
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_vary%s_nearest',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_vary%s_nearest',opt.fNamePostFix, opt.addOffsetParam, opt.vary.perturbOrigPRFs)),[],[1 300],'.',1);
    
end


%% Plot time series
if ~opt.vary.perturbOrigPRFs
    load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'megStimConditions.mat'), 'triggers');
    
    % Define time scale
    [~, nSensors,nEpochs] = size(allData);
    epochLengthSeconds = diff(opt.meg.epochStartEnd);
    t = (0:nEpochs-1) .* epochLengthSeconds;
    
    % Define blink and blank blocks
    blinkIdx = (triggers.stimConditions(1:length(t))==20);
    blankIdx = (triggers.stimConditions(1:length(t))==10);
    blink_t = t(blinkIdx);
    blank_t = t(blankIdx);
    
    averageAllData = nanmean(allData,1);
    averageAllPredictions = nanmean(allPredictions,1);
    
    % Plot individual sensor time series
    fH3 = figure(3); set(fH3, 'Position', [1000,420,1269,918]);
    for s = 1:nSensors
        
        % Compute y limits
        tmp_yl = max(abs([min(averageAllData(s,:)), max(averageAllData(s,:))])).*10^14;
        if (tmp_yl > 3)
            yl = [-1*tmp_yl, tmp_yl].*10^-14;
        else
            yl = [-3,3].*10^-14;
        end
        ylim(yl); xlim([0, max(t)])
        
        clf,
        averageAllData(s,blinkIdx) = NaN;
        
        % Plot the blank and blink periods
        for bt = 1:length(blink_t)
            patch([blink_t(bt),blink_t(bt)+epochLengthSeconds blink_t(bt)+epochLengthSeconds blink_t(bt)],[yl(1),yl(1),yl(2),yl(2)],[0.5 0.5 0.5],'FaceAlpha', 0.2, 'LineStyle','none');
        end
        
        for bt = 1:length(blank_t)
            patch([blank_t(bt),blank_t(bt)+epochLengthSeconds blank_t(bt)+epochLengthSeconds blank_t(bt)],[yl(1),yl(1),yl(2),yl(2)],[0.5 0.5 0.5],'FaceAlpha', 0.7, 'LineStyle','none');
        end
        hold on;
        plot(t, zeros(size(timepoints)), 'k');
        plot(t, averageAllData(s,:), 'ko:', 'LineWidth',3);
        plot(t, averageAllPredictions(s,:), '-', 'Color',[1 0.45 0.45], 'LineWidth',3);
        ylim(yl); xlim([0, max(t)])
        title(sprintf('Group average sensor %d, var expl %1.2f', s, groupVarExpl(s)));
        set(gca, 'TickDir', 'out', 'FontSize', 20, 'TickLength',[0.010 0.010],'LineWidth',3); box off;
        xlabel('Time (s)'); ylabel('Phase-referenced 10 Hz SSVEF response (fT)');
        
        l = findobj(gca, 'Type','Line');
        legend(l([2,1]), {'Observed', 'Predicted'}, 'Location', 'NorthEastOutside', 'FontSize', 25); legend boxoff;
        
        figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_sensor_%d_%s_Offset%d',s, opt.fNamePostFix, opt.addOffsetParam)),[],0,'.',1);
        figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_sensor_%d_%s_Offset%d',s, opt.fNamePostFix, opt.addOffsetParam)),[],[1 300],'.',1);
        
    end
    
    %     % Plot individual subjects' time series on top of eachother + mean
    %     cmap = lines(length(subjIDs));
    %     fH4 = figure(4); set(fH4, 'Position', [1000,420,1269,918]);
    %     for s = 1:nSensors
    %
    %         subplot(211); cla;
    %         plot(t, zeros(size(t)), 'k'); hold on;
    %
    %         subplot(212); cla;
    %         plot(t, zeros(size(t)), 'k'); hold on;
    %
    %         for subj = 1:length(subjIDs)
    %             subplot(211)
    %             plot(t, squeeze(allPredictions(subj,:,s)), '-', 'Color', cmap(subj,:), 'LineWidth',1);
    %
    %             subplot(212)
    %             plot(t, squeeze(allData(subj,:,s)), '-', 'Color', cmap(subj,:), 'LineWidth',1);
    %         end
    %
    %         subplot(211);
    %         plot(t, groupAvePredictionScaled(:,s).*10^-14, 'r', 'LineWidth',4);
    %         title(sprintf('Individual predictions sensor %d', s));
    %         set(gca, 'TickDir', 'out', 'FontSize', 20, 'TickLength',[0.010 0.010]); box off;
    %         xlabel('Time (s)'); ylabel('Phase-ref 10 Hz SSVEF (fT)');
    %
    %         subplot(212)
    %         plot(t, groupAveData(:,s), 'k', 'LineWidth',4);
    %         title(sprintf('Individual data sensor %d', s));
    %         set(gca, 'TickDir', 'out', 'FontSize', 20, 'TickLength',[0.010 0.010]); box off;
    %         xlabel('Time (s)'); ylabel('Phase-ref 10 Hz SSVEF (fT)');
    %
    %
    %         figurewrite(fullfile(saveDir, sprintf('IndivSubjectsData_sensor_%d_%s_Offset%d',s, opt.fNamePostFix,opt.addOffsetParam)),[],0,'.',1);
    %         figurewrite(fullfile(saveDir, sprintf('IndivSubjectsData_sensor_%d_%s_Offset%d',s, opt.fNamePostFix, opt.addOffsetParam)),[],[1 300],'.',1);
    %
    %     end
    
end