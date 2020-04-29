function [] = mprf_CompareGroupAverageDataVsPrediction(dirPth, opt, whichFigure, sensorsToAverage)

subjIDs = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058', ...
    'wlsubj068', 'wlsubj070', 'wlsubj081', 'wlsubj106', ...
    'wlsubj109', 'wlsubj111'};

saveSubDir = sprintf('figure%d_%s',whichFigure, opt.subfolder);
saveDir = fullfile(dirPth.finalFig.savePthAverage,saveSubDir, 'GroupAvePrediction');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

%% Preallocate space
nEpochs  = 140;
nSensors = length(opt.meg.dataChan);

if ~opt.vary.perturbOrigPRFs
    nrVariations = 1;
elseif strcmp(opt.vary.perturbOrigPRFs, 'position')
    nrVariations = length(opt.vary.position);
elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
    nrVariations = length(opt.vary.size);
end

allPredictions = NaN(length(subjIDs), nEpochs, nSensors, nrVariations);
allData        = NaN(length(subjIDs), nEpochs, nSensors, nrVariations);

%% Concatenate data
for subj = 1:length(subjIDs)
    
    subjectID = subjIDs{subj};
    
    dirPth = loadPaths(subjectID);
    
    % Load sensor predictions
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp', 'predMEGResponseScaled'));
    
    % Load phase-referenced SSVEF data
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp', 'meanPhRefAmp10Hz'));
    
    allPredictions(subj,:,:,:) = predMEGResponseScaled;
    allData(subj,:,:,:) = meanPhRefAmp10Hz;
end

%% Predict group average data from group average predictions

% Take the average across subjects
groupAvePrediction = squeeze(nanmean(allPredictions,1));
groupAveData       = squeeze(nanmean(allData,1));

% Preallocate space
groupAvePredictionScaled = NaN(nEpochs, nSensors, nrVariations);
groupVarExpl             = NaN(nSensors, nrVariations);

% Define timepoints
timepoints = 1:size(allData,2);
femtoScaleFactor = 10^14;

for v = 1:nrVariations
    
    for s = 1:nSensors
        
        % Identify and remove nans
        meanNanMask = isnan(groupAveData(:,s,v));
        thisGroupAvePredictionMasked = groupAvePrediction(~meanNanMask,s,v) .*femtoScaleFactor;
        thisGroupAveDataMasked       = groupAveData(~meanNanMask,s,v) .*femtoScaleFactor;
        
        % Create predictions
        if opt.addOffsetParam
            
            groupAveX = [ones(size(thisGroupAvePredictionMasked,1),1), thisGroupAvePredictionMasked];
            
            % Regress out predictions
            groupFitBeta = groupAveX \ thisGroupAveDataMasked;
            
            % Compute scaled predictions with betas
            groupAvePredictionScaled(~meanNanMask,s,v) =  thisGroupAvePredictionMasked * groupFitBeta(2) + groupFitBeta(1);
        else
            groupAveX = thisGroupAvePredictionMasked;
            
            % Regress out predictions
            groupFitBeta = groupAveX \ thisGroupAveDataMasked;
            
            % Compute scaled predictions with betas
            groupAvePredictionScaled(~meanNanMask,s,v) =  groupAveX * groupFitBeta';
        end
        
        % Compute coefficient of determination:
        groupVarExpl(s,v) = computeCoD(thisGroupAveDataMasked,groupAvePredictionScaled(~meanNanMask,s,v));
        
    end
    
end

% Remove dummy dimension
groupVarExpl             = squeeze(groupVarExpl);
groupAvePredictionScaled = squeeze(groupAvePredictionScaled);


%% PLOTTING
if ~opt.vary.perturbOrigPRFs
    
    % Plot var explained map
    clims = [0 max(groupVarExpl)];
    fH1 = figure; clf;
    megPlotMap(groupVarExpl,clims,fH1, 'parula', ...
        'Average Group Fit Variance Explained', [],[], 'interpMethod', 'v4');
    c = colorbar;
    c.Location = 'eastoutside';
    c.Box = 'off';
    c.TickDirection = 'out';
    c.TickLength = [0.010 0.010];
    c.FontSize = 12;
    
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d',opt.fNamePostFix, opt.addOffsetParam)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d',opt.fNamePostFix, opt.addOffsetParam)),[],[1 300],'.',1);
    
    fH1 = figure; clf;
    megPlotMap(groupVarExpl,clims,fH1, 'parula', ...
        'Average Group Fit Variance Explained', [],[], 'interpmethod', 'nearest');
    c = colorbar;
    c.Location = 'eastoutside';
    c.Box = 'off';
    c.TickDirection = 'out';
    c.TickLength = [0.010 0.010];
    c.FontSize = 12;
    
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_nearest',opt.fNamePostFix, opt.addOffsetParam)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_Offset%d_nearest',opt.fNamePostFix, opt.addOffsetParam)),[],[1 300],'.',1);
    
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
    sensorLoc = selectSensorsToAverage(opt, dirPth, saveDir, groupVarExpl', sensorsToAverage);
    
    % Select data to plot
    if strcmp(sensorsToAverage, 'top10Positive')
        curSubjLocs = sensorLoc;
        for ii = 1:size(curSubjLocs,1)
            curSensorLoc = curSubjLocs(ii,:);
            dataToPlot(ii,1:sum(~isnan(curSensorLoc))) = squeeze(groupVarExpl(curSensorLoc(~isnan(curSensorLoc)),ii));
        end
    else
        dataToPlot = squeeze(groupVarExpl(sensorLoc,:))';
    end
    
    
    figure(99); clf;
    meandataToPlot = 100.*nanmean(dataToPlot, 2);
    % Plot mean with shaded error bar using 'patch' function
    se   = 100.*nanstd(dataToPlot, [],2) ./ sqrt(size(dataToPlot,2));
    ci   = 1 .* se; % use zsore=1 for 68% CI or zcore=1.95 for 95%
    lo = meandataToPlot - ci;
    hi = meandataToPlot + ci;
    
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
        
        figure(fH1); subplot(2,ceil(nrVariations/2),ii)
        megPlotMap(groupVarExpl(:,ii)',clims,fH1, 'parula', ...
            subplotLabels{ii}, [],[], 'interpMethod', 'v4');
        c = colorbar;
        c.Location = 'eastoutside';
        c.Box = 'off';
        c.TickDirection = 'out';
        c.TickLength = [0.010 0.010];
        c.FontSize = 12;
        pos = c.Position; set(c, 'Position', [pos(1)+0.03 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
        
        figure(fH2); subplot(3,ceil(nrVariations/3),ii)
        megPlotMap(groupVarExpl(:,ii)',clims,fH2, 'parula', ...
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
    [nEpochs, nSensors] = size(meanPhRefAmp10Hz);
    epochLengthSeconds = diff(opt.meg.epochStartEnd);
    t = (0:nEpochs-1) .* epochLengthSeconds;
    
    % Define blink and blank blocks
    blinkIdx = (triggers.stimConditions(1:length(t))==20);
    blankIdx = (triggers.stimConditions(1:length(t))==10);
    blink_t = t(blinkIdx);
    blank_t = t(blankIdx);
    
    
    % Plot individual sensor time series
    fH3 = figure(3); set(fH3, 'Position', [1000,420,1269,918]);
    for s = 1:nSensors
        
        % Compute y limits
        tmp_yl = max(abs([min(groupAveData(:,s)), max(groupAveData(:,s))])).*10^14;
        if (tmp_yl > 3)
            yl = [-1*tmp_yl, tmp_yl].*10^-14;
        else
            yl = [-3,3].*10^-14;
        end
        ylim(yl); xlim([0, max(t)])
        
        clf,
        groupAveData(blinkIdx,s) = NaN;
        
        % Plot the blank and blink periods
        for bt = 1:length(blink_t)
            patch([blink_t(bt),blink_t(bt)+epochLengthSeconds blink_t(bt)+epochLengthSeconds blink_t(bt)],[yl(1),yl(1),yl(2),yl(2)],[0.5 0.5 0.5],'FaceAlpha', 0.2, 'LineStyle','none');
        end
        
        for bt = 1:length(blank_t)
            patch([blank_t(bt),blank_t(bt)+epochLengthSeconds blank_t(bt)+epochLengthSeconds blank_t(bt)],[yl(1),yl(1),yl(2),yl(2)],[0.5 0.5 0.5],'FaceAlpha', 0.7, 'LineStyle','none');
        end
        hold on;
        plot(t, zeros(size(timepoints)), 'k');
        plot(t, groupAveData(:,s), 'ko:', 'LineWidth',3);
        plot(t, groupAvePredictionScaled(:,s).*10^-14, '-', 'Color',[1 0.45 0.45], 'LineWidth',3);
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
