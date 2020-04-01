function [] = mprf_CompareGroupAverageDataVsPrediction(dirPth, opt)

subjIDs = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058', ...
        'wlsubj068', 'wlsubj070', 'wlsubj081', 'wlsubj106', ...
        'wlsubj109', 'wlsubj111'};

saveSubDir = ['figure1_' opt.regressionType];
saveDir = fullfile(dirPth.finalFig.savePthAverage,saveSubDir, 'GroupAvePrediction');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

nEpochs = 140;
nSensors = 157;


%% Concatenate data
allPredictions = NaN(length(subjIDs), nEpochs, nSensors);
allData        = NaN(length(subjIDs), nEpochs, nSensors);

for subj = 1:length(subjIDs)
    
    subjectID = subjIDs{subj};
    
    dirPth = loadPaths(subjectID);
    
    % Load sensor predictions
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp', 'predMEGResponseScaled'));
    
    % Load phase-referenced SSVEF data
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp', 'meanPhRefAmp10Hz'));
    
    allPredictions(subj,:,:) = predMEGResponseScaled;
    allData(subj,:,:) = meanPhRefAmp10Hz;
end


%% Predict group average data from group average predictions
groupAvePrediction = squeeze(nanmean(allPredictions,1));
groupAveData       = squeeze(nanmean(allData,1));


groupAveDataScaled = NaN(size(groupAvePrediction));
groupAvePredictionScaled = NaN(size(groupAvePrediction));
timepoints = 1:size(allData,2);

for s = 1:nSensors
    
    % Identify and remove nans
    meanNanMask = isnan(groupAveData(:,s));
    thisGroupAvePredictionMasked = groupAvePrediction(~meanNanMask,s) .*10^14;
    thisGroupAveDataMasked       = groupAveData(~meanNanMask,s) .*10^14;
    
    % Create predictions
    if strcmp(opt.regressionType, 'NoOffset')
        meanX = thisGroupAvePredictionMasked;
        
        % Regress out predictions
        meanB = meanX \ thisGroupAveDataMasked;
        
        % Compute scaled predictions with betas
        groupAvePredictionScaled(~meanNanMask,s) =  meanX * meanB';
        
    else
        meanX = [ones(size(thisGroupAvePredictionMasked,1),1), thisGroupAvePredictionMasked];
    
        % Regress out predictions
        meanB = meanX \ thisGroupAveDataMasked;
    
        % Compute scaled predictions with betas
        groupAvePredictionScaled(~meanNanMask,s) =  thisGroupAvePredictionMasked * meanB(2) + meanB(1);
    end
    
    % Compute coefficient of determination:
    groupVarExpl(s) = 1 - (sum( (thisGroupAveDataMasked - groupAvePredictionScaled(~meanNanMask,s)).^2) ...
        ./ sum(thisGroupAveDataMasked.^2));
    
end


% Plot var explained map
clims = [0 max(groupVarExpl)];
fH1 = figure; clf;
megPlotMap(groupVarExpl,clims,fH1, 'parula', ...
    'Average Group Fit Variance Explained', [],[]);
c = colorbar;
c.Location = 'eastoutside';
c.Box = 'off';
c.TickDirection = 'out';
c.TickLength = [0.010 0.010];
c.FontSize = 12;

figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_%s',opt.fNamePostFix, opt.regressionType)),[],0,'.',1);
figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_%s',opt.fNamePostFix, opt.regressionType)),[],[1 300],'.',1);


%% Plot time series

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
fH2 = figure(2); set(fH2, 'Position', [1000,420,1269,918]);
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
    
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_sensor_%d_%s_%s',s, opt.fNamePostFix, opt.regressionType)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_sensor_%d_%s_%s',s, opt.fNamePostFix, opt.regressionType)),[],[1 300],'.',1);
    
end

% Plot individual subjects' time series on top of eachother + mean
cmap = lines(length(subjIDs));
fH3 = figure(3); set(fH3, 'Position', [1000,420,1269,918]);
for s = 1:nSensors
    
    subplot(211); cla;
    plot(t, zeros(size(t)), 'k'); hold on;
    
    subplot(212); cla;
    plot(t, zeros(size(t)), 'k'); hold on;
    
    for subj = 1:length(subjIDs)
        subplot(211)
        plot(t, squeeze(allPredictions(subj,:,s)), '-', 'Color', cmap(subj,:), 'LineWidth',1);
        
        subplot(212)
        plot(t, squeeze(allData(subj,:,s)), '-', 'Color', cmap(subj,:), 'LineWidth',1);
    end
    
    subplot(211);
    plot(t, groupAvePredictionScaled(:,s).*10^-14, 'r', 'LineWidth',4);
    title(sprintf('Individual predictions sensor %d', s));
    set(gca, 'TickDir', 'out', 'FontSize', 20, 'TickLength',[0.010 0.010]); box off;
    xlabel('Time (s)'); ylabel('Phase-ref 10 Hz SSVEF (fT)');
    
    subplot(212)
    plot(t, groupAveData(:,s), 'k', 'LineWidth',4);
    title(sprintf('Individual data sensor %d', s)); 
    set(gca, 'TickDir', 'out', 'FontSize', 20, 'TickLength',[0.010 0.010]); box off;
    xlabel('Time (s)'); ylabel('Phase-ref 10 Hz SSVEF (fT)');
      
    
    figurewrite(fullfile(saveDir, sprintf('IndivSubjectsData_sensor_%d_%s_%s',s, opt.fNamePostFix,opt.regressionType)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('IndivSubjectsData_sensor_%d_%s_%s',s, opt.fNamePostFix, opt.regressionType)),[],[1 300],'.',1);
    
end
