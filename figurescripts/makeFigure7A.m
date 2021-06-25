function makeFigure7A(sensorsToAverage, summaryMetric, opt)
% Function to make sensorwise and bootstrapped sensorwise average across subjects,
% plotting variance explained by the model as a function of polar angle
% rotations around the fovea of the original estimated pRF centers.

% Check inputs
if ~exist('sensorsToAverage', 'var') || isempty(sensorsToAverage)
    sensorsToAverage = 'top10';
end

if ~exist('summaryMetric', 'var') || isempty(summaryMetric)
    summaryMetric = 'meanVE';
end

if ~exist('opt', 'var') || isempty(opt)
    opt = getOpts('perturbOrigPRFs','position');
end

doDemeanVE = false;

% Define all subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

% Define the range of rotations
range   = opt.vary.position;

% Allocate space
varexpl = NaN(length(subjects),length(range), 157);
sensorLoc = cell(length(subjects),1);

% Set up figures 
fH1   = figure(1); clf; set(fH1, 'Color', 'w', 'Position', [1, 592, 838, 746]);

% Plotting params
lineColorSat = repmat(linspace(0.3,0.9,10), [3 1]);

for s = 1:length(subjects)
    
    % Get subject name and directories
    subjectID = subjects{s};
    dirPth = loadPaths(subjectID);
    
    [pth, ~] = fileparts(dirPth.model.saveFigPth);
    saveSubDir = ['SupplementaryFigure9_varyPosition'];
    saveDir = fullfile(pth, 'finalfig', saveSubDir);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    
    % Load variance explained file
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp', 'meanVarExpl'));
    varexpl(s,:,:) = meanVarExpl;
    
    % What sensors are we averaging?
    sensorLoc{s} = selectSensorsToAverage(opt, dirPth, saveDir, squeeze(varexpl(s,:,:)), sensorsToAverage);
    
    % Select those data corresponding to sensors
    if strcmp(sensorsToAverage, 'top10Positive')
        curSubjLocs = sensorLoc{s};
        for ii = 1:size(curSubjLocs,1)
            curSensorLoc = curSubjLocs(ii,:);
            thisSubjectSensorData(ii,1:sum(~isnan(curSensorLoc))) = squeeze(varexpl(s,ii,curSensorLoc(~isnan(curSensorLoc))));
        end
    else
        thisSubjectSensorData = squeeze(varexpl(s,:,sensorLoc{s}));
    end
    
    % Compute summary metrics of variance explained across selected sensors
    meanSelectedSensors(s,:) = 100*nanmean(thisSubjectSensorData,2);
    
    if strcmp(summaryMetric, 'meanVE')
        dataToPlot = meanSelectedSensors;
        yLabel = 'Variance explained (%)';
        yl = [-10 50];
        % Rescale if subjects fall outside ylimit
        if max(dataToPlot(:))>yl(2)
            yl = [yl(1) max(dataToPlot(:))+10];
        end
        if min(dataToPlot(:))<yl(1)
            yl = [min(dataToPlot(:))-10 yl(2)];
        end
    elseif strcmp(summaryMetric, 'percentChangeVE')
        dataToPlot(s,:) = 100*((meanSelectedSensors(s,:) - mean(meanSelectedSensors(s,:)))./mean(meanSelectedSensors(s,:)));
        yl = [-100 100];
        yLabel = 'Percent change variance explained (%)';
    elseif strcmp(summaryMetric, 'zscoreVE')
        dataToPlot(s,:) = zscore(meanSelectedSensors(s,:));
        yl = [-3 3];
        yLabel = 'Z-scored variance explained (%)';
    end
    
    figure(fH1);
    plot(range,dataToPlot(s,:),'Color', lineColorSat(:,s), 'Linewidth',2); hold on;
    
end

%% Plot average across subjects on top of figure with individual lines
if doDemeanVE
    mn = nanmean(dataToPlot,2);
    averageDataToPlot = nanmean(dataToPlot-mn,1) + mean(mn,1);
else
    averageDataToPlot = nanmean(dataToPlot,1);
end



figure(fH1);
plot(range,averageDataToPlot,'r','Linewidth',5); hold on;
plot(range, zeros(size(averageDataToPlot)), 'k', 'LineWidth', 1);
plot([0 0], [min(yl), max(yl)], 'k');

yl = [0 50];
% Add labels and make pretty
set(gca,'TickDir', 'out');
xlabel('Rotation angle from original pRF position (deg)');
set(gca,'XTick', range,'XTickLabel',rad2deg(range), 'YLim', yl, 'XLim', [range(1),range(end)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
title('sensorwise Average variance explained by modelfit: Vary Position');
ylabel(yLabel);
box off;

%% Figure 2 with bootstrapping data across subjects

colorCIPatch = [0.7 0.7 0.7];

% Bootstrap with 10,000 iterations
nBoot = 10000;
if doDemeanVE
    mnSubs = mean(dataToPlot,2,'omitnan');
    BootStrappedData = bootstrp(nBoot, @(x) mprf_averageVar(x,dataToPlot-mnSubs), (1:size(dataToPlot,1)));
else
    BootStrappedData = bootstrp(nBoot, @(x) mprf_averageVar(x,dataToPlot), (1:size(dataToPlot,1)));
end

pct1 = 100 * (0.32/2);
pct2 = 100 - pct1;
lo = prctile(BootStrappedData,pct1); % 16th percentile
hi = prctile(BootStrappedData,pct2); % 84th percentile

if doDemeanVE
    aveBoot = nanmean(BootStrappedData,1);
    aveBoot = aveBoot + mnSubs;
    lo = lo+mnSubs;
    hi = hi+mnSubs;
else
    % Average of bootstrapped data
    aveBoot = nanmean(BootStrappedData,1);
end

% Get p value from bootstrapped data
% origIdx = find(range==0);
% [~,peakIdx] = max(BootStrappedData(:,1:origIdx+2),[],2);
% peakBelowOrigPRF = sum(peakIdx~=origIdx);
% p = 2*(.5-abs(.5- (peakBelowOrigPRF/nBoot)));
% fprintf('Mean different from original pRF position of %dx bootstrapped p value: %1.3f \n',nBoot, p)

% Plot bootstrapped data!
yl = [5 30];
fH2 = figure(2); clf; set(gcf,'position',[66,1,1855,1001]);
plot(range, zeros(size(aveBoot)), 'k', 'LineWidth', 1); hold on;
patch([range, fliplr(range)], [lo, fliplr(hi)],colorCIPatch, 'FaceAlpha', 0.5, 'LineStyle',':');
plot(range,aveBoot,'r','Linewidth',5);
plot([0 0], [min(yl), max(yl)], 'k');

% Add labels and make pretty
set(gca,'TickDir', 'out');
xlabel('Rotation angle from original pRF position (deg)');
set(gca,'XTick', range,'XTickLabel',rad2deg(range), 'YLim', yl, 'XLim', [range(1),range(end)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
title('Fit-then-average variance explained by modelfit: Vary Position');
ylabel(yLabel);
box off;

if opt.saveFig
   
    saveSubDir = ['Figure7A_fitThenAverage_' opt.subfolder];
    saveDir = fullfile(dirPth.finalFig.savePthAverage, saveSubDir, sensorsToAverage);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    fprintf('\n(%s): Saving Figure 7A in %s\n',mfilename, saveDir);
    print(fH2, fullfile(saveDir, sprintf('SupplFigure7A_%s_varyPositionFitThenAverageSummary%s_%s_demean%d', dirPth.subjID, opt.fNamePostFix, sensorsToAverage, doDemeanVE)), '-dpdf');
    figurewrite(fullfile(saveDir, sprintf('SupplFigure7A_%s_varyPositionFitThenAverageSummary%s_%s_demean%d', dirPth.subjID, opt.fNamePostFix, sensorsToAverage,doDemeanVE)), [],[1 300],'.',1);

end

return