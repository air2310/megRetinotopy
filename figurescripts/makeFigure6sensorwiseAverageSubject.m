function makeFigure6sensorwiseAverageSubject(sensorsToAverage, summaryMetric, opt)
% Function to make sensorwise average and bootstrapped average across subjects,
% plotting variance explained by the model as a function of scale factor
% of the original estimated pRF size.

% Check inputs
if ~exist('sensorsToAverage', 'var') || isempty(sensorsToAverage)
    sensorsToAverage = 'top10';
end

if ~exist('summaryMetric', 'var') || isempty(summaryMetric)
    summaryMetric = 'meanVE';
end

if ~exist('opt', 'var') || isempty(opt)
    opt = getOpts;
end


% Define subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

% Define the range of rotations
range   = opt.vary.size;

% Allocate space
varexpl = NaN(length(subjects),length(range), 157);
sensorLoc = cell(length(subjects),1);

% Get sensor locations in the back
load(which('meg160_example_hdr.mat'))
layout = ft_prepare_layout([],hdr);
xpos = layout.pos(1:157,1);
ypos = layout.pos(1:157,2);

% Figure specs
fH1   = figure(1); clf; set(fH1, 'Color', 'w', 'Position', [1, 592, 838, 746]);
lineColorSat = repmat(linspace(0.3,0.9,10), [3 1]);

for s = 1:length(subjects)
    % Get subject name and directories
    subjectID = subjects{s};
    dirPth = loadPaths(subjectID);
    
    % Load variance explained file
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp', 'meanVarExpl'));
    varexpl(s,:,:) = meanVarExpl;
    
    % What sensors are we averaging?
    sensorLoc{s} = selectSensorsToAverage(opt, dirPth, saveDir, squeeze(varexpl(s,:,:)), sensorsToAverage);
    
    % Select data for sensors over the back of the head
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
        yl = [-10 50];
        % Rescale if subjects fall outside ylimit
        if max(dataToPlot(:))>yl(2)
            yl = [yl(1) max(dataToPlot(:))+10];
        end
        if min(dataToPlot(:))<yl(1)
            yl = [min(dataToPlot(:))-10 yl(2)];
        end
        yLabel = 'Variance explained (%)';
    elseif strcmp(summaryMetric, 'percentChangeVE')
        dataToPlot(s,:) = 100*((meanSelectedSensors(s,:) - mean(meanSelectedSensors(s,:)))./mean(meanSelectedSensors(s,:)));
        yl = [-100 100];
        yLabel = 'Percent change variance explained (%)';
    elseif strcmp(summaryMetric, 'zscoreVE')
        dataToPlot(s,:) = zscore(meanSelectedSensors(s,:));
        yl = [-3 3];
        yLabel = 'Z-scored variance explained (%)';
    end
    
    % Compute mean and standard error of variance explained across selected sensors
    figure(fH1);
    plot(range,meanSelectedSensors(s,:),'Color', lineColorSat(:,s), 'Linewidth',2); hold on;

end % subjects

%% Plot average across subjects on top of figure with individual lines

% Average across subjects
averageDataToPlot = nanmean(dataToPlot,1);

figure(fH1);
plot(range,averageDataToPlot,'r','Linewidth',6);
plot([1 1], [min(yl), max(yl)], 'k');

% Add labels and make pretty
set(gca,'TickDir', 'out');
xlabel('Scale factor of original pRF size');
set(gca,'XTick', range,'XTickLabel',range, 'YLim', yl, 'XLim', [range(1),range(end)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20, 'XScale', 'log'); axis square;
title('Sensorwise average of variance explained by modelfit: Vary Size');
ylabel(yLabel);
box off;

%% Figure 3 with bootstrapping data across subjects
% Bootstrap average variance explained across subjects with 10,000 iterations
nBoot = 10000;
BootStrappedData = bootstrp(nBoot, @(x) mprf_averageVar(x,dataToPlot), (1:size(dataToPlot,1)));
pct1 = 100 * (0.32/2);
pct2 = 100 - pct1;
lo = prctile(BootStrappedData,pct1); % 16th percentile
hi = prctile(BootStrappedData,pct2); % 84th percentile

% Get p value from bootstrapped data
origIdx = find(range==1);
[~,peakIdx] = max(BootStrappedData(:,1:origIdx+2),[],2);
peakBelowOrigPRF = sum(peakIdx<origIdx);
p = 2*(.5-abs(.5- (peakBelowOrigPRF/nBoot)));
fprintf('Mean below original pRF size of %dx bootstrapped p value: %1.3f \n',nBoot, p)

% Average of bootstrapped data
aveBoot = nanmean(BootStrappedData,1);

fH2 = figure(2); clf; set(gcf,'Position',[1, 592, 838, 746]);
plot(range,zeros(size(aveBoot)),'k','Linewidth',1); hold on;
patch([range, fliplr(range)], [lo, fliplr(hi)],colorCIPatch, 'FaceAlpha', 0.5, 'LineStyle',':');
plot(range,aveBoot,'r','Linewidth',5);
plot([1 1], [min(lo), max(hi)], 'k');

% Add labels and make pretty
set(gca,'TickDir', 'out');
xlabel('Scale factor of original pRF size');
set(gca,'XTick', range,'XTickLabel',range, 'YLim', [min(lo), max(hi)], 'XLim', [range(1),range(end)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20, 'XScale', 'log'); axis square;
title('Bootstrapped average variance explained by modelfit: Vary Size');
ylabel(yLabel);
box off;


return