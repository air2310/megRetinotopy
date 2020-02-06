% mprf_CompareGroupAverageDataVsPrediction

opt = getOpts('saveFig',1,'verbose',0, 'fullSizeMesh', 1, 'perturbOrigPRFs',0); % see getOpts function for more options

subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

nEpochs = 140;
nSensors = 157;

allPredictions = NaN(length(subjects), nEpochs, nSensors);
allData        = NaN(length(subjects), nEpochs, nSensors);


for subj = 1:length(subjects)
    
    subjectID = subjects{subj};
    
    dirPth = loadPaths(subjectID);
    
    % Load sensor predictions
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp', 'meanPredResponse'));
    
    % Load beta weights and ref phases
    %     load(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp', 'bestBetas_1'));
    
    % Load phase-referenced SSVEF data
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp', 'meanPhRefAmp10Hz'));
    
    allPredictions(subj,:,:) = meanPredResponse;
    allData(subj,:,:) = meanPhRefAmp10Hz;
end




groupAvePrediction = squeeze(nanmean(allPredictions,1));
groupAveData       = squeeze(nanmean(allData,1));


groupAveDataScaled = NaN(size(groupAvePrediction));
t = 1:size(allData,2);


saveSubDir = 'figureGroupAvePrediction';
saveDir = fullfile(dirPth.finalFig.savePthAverage,'figureSuperSubject',saveSubDir);
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end



for s = 1:nSensors
    
    % Identify and remove nans
    meanNanMask = isnan(groupAvePrediction(:,s));
    thisGroupAvePredictionMasked = groupAvePrediction(~meanNanMask,s);
    thisGroupAveDataMasked       = groupAveData(~meanNanMask,s);
    
    % Create predictions
%     meanX = [ones(size(thisGroupAvePredictionMasked)) thisGroupAvePredictionMasked];
    meanX = thisGroupAvePredictionMasked;
    
    % Regress out predictions
    meanB = meanX \ thisGroupAveDataMasked;
    
    % Compute scaled predictions with betas
    groupAvePredictionScaled(~meanNanMask,s) =  meanX * meanB';
    
    % Compute coefficient of determination:
    groupVarExpl(s) = 1 - (var(thisGroupAveDataMasked - groupAvePredictionScaled(~meanNanMask,s)) ...
        ./ var(thisGroupAveDataMasked));
    
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

figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_zeroblankpred',opt.fNamePostFix)),[],0,'.',1);
figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_VarExplMesh_%s_zeroblankpred',opt.fNamePostFix)),[],[1 300],'.',1);


% Plot individual sensor time series
fH2 = figure(2); set(fH2, 'Position', [1000,420,1269,918]);
for s = 1:nSensors
    clf,
    plot(t, zeros(size(t)), 'k'); hold on;
    plot(t, groupAveData(:,s), 'ko:', 'LineWidth',2);
    plot(t, groupAvePredictionScaled(:,s), 'r-', 'LineWidth',2);
    title(sprintf('Group average sensor %d, var expl %1.2f', s, groupVarExpl(s)));
    set(gca, 'TickDir', 'out', 'FontSize', 20, 'TickLength',[0.010 0.010]); box off;
    xlabel('Time (s)'); ylabel('Phase-referenced 10 Hz SSVEF response (fT)');
    
    % Compute y limits
    tmp_yl = max(abs([min(groupAveData(:,s)), max(groupAveData(:,s))])).*10^14;
    if (tmp_yl > 3)
        yl = [-1*tmp_yl, tmp_yl].*10^-14;
    else
        yl = [-3,3].*10^-14;
    end
    ylim(yl); xlim([0, max(t)])
    
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_sensor_%d_%s_zeroblankpred',s, opt.fNamePostFix)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('GroupAverageFit_sensor_%d_%s_zeroblankpred',s, opt.fNamePostFix)),[],[1 300],'.',1);
    
end

% Plot individual subjects' time series on top of eachother + mean
cmap = lines(length(subjects));
figure(2); clf;
for s = 1:nSensors
    
    subplot(211); cla;
    plot(t, zeros(size(t)), 'k'); hold on;
    
    subplot(212); cla;
    plot(t, zeros(size(t)), 'k'); hold on;
    
    for subj = 1:length(subjects)
        subplot(211)
        plot(t, squeeze(allPredictions(subj,:,s)), '-', 'Color', cmap(subj,:), 'LineWidth',1);
        
        subplot(212)
        plot(t, squeeze(allData(subj,:,s)), '-', 'Color', cmap(subj,:), 'LineWidth',1);
    end
    
    subplot(211);
    plot(t, groupAvePrediction(:,s), 'r', 'LineWidth',4);
    title(sprintf('Individual predictions sensor %d', s));
    set(gca, 'TickDir', 'out', 'FontSize', 20, 'TickLength',[0.010 0.010]); box off;
    xlabel('Time (s)'); ylabel('Phase-ref 10 Hz SSVEF (fT)');
    
    subplot(212)
    plot(t, groupAveData(:,s), 'k', 'LineWidth',4);
    title(sprintf('Individual data sensor %d', s)); 
    set(gca, 'TickDir', 'out', 'FontSize', 20, 'TickLength',[0.010 0.010]); box off;
    xlabel('Time (s)'); ylabel('Phase-ref 10 Hz SSVEF (fT)');
      
    
    figurewrite(fullfile(saveDir, sprintf('IndivSubjectsData_sensor_%d_%s_zeroblankpred',s, opt.fNamePostFix)),[],0,'.',1);
    figurewrite(fullfile(saveDir, sprintf('IndivSubjectsData_sensor_%d_%s_zeroblankpred',s, opt.fNamePostFix)),[],[1 300],'.',1);
    
end
