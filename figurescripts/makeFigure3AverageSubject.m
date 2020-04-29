function makeFigure3AverageSubject(dirPth, opt, sensorsToAverage, summaryMetric)
% Function to make average across subjects for Figure 2 from manuscript,
% plotting variance explained by the model as a function of polar angle
% rotations around the fovea of the original estimated pRF centers.

% Check inputs
if ~exist('summaryMetric', 'var') || isempty(summaryMetric)
    summaryMetric = 'meanVE';
end

% Define folder to save figures
if opt.saveFig
    saveSubDir = ['figure3_' opt.subfolder];
    saveDir = fullfile(dirPth.finalFig.savePthAverage, saveSubDir);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
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

% Figure for all subjects plotted separately
if  strcmp(summaryMetric, 'meanVE')
    sz = get(0, 'screensize');
    fH2 = figure; set(gcf,'Position', sz);
    colorCIPatch = [0.5 0.5 0.5];
end

for s = 1:length(subjects)
    
    subjectID = subjects{s};
    
    dirPth = loadPaths(subjectID);
    
    % Load variance explained file
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp', 'meanVarExpl'));
    
    varexpl(s,:,:) = meanVarExpl;
    
    % What sensors are we averaging?
    sensorLoc{s} = selectSensorsToAverage(opt, dirPth, saveDir, squeeze(varexpl(s,:,:)), sensorsToAverage);
    
    % Plot data for sensors over the back of the head
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
            yl = [0 60];
            % Rescale if subjects fall outside ylimit
            if max(dataToPlot(:))>yl(2)
                yl = [yl(1) max(dataToPlot(:))+10];
            end
            if min(dataToPlot(:))<yl(1)
                yl = [-50 yl(2)];
            end
            yLabel = 'Variance explained (%)';
        elseif strcmp(summaryMetric, 'percentChangeVE')
            percentdiff(s,:) = 100*((meanSelectedSensors(s,:) - mean(meanSelectedSensors(s,:)))./mean(meanSelectedSensors(s,:)));
            dataToPlot = percentdiff;
            yl = [-100 100];
            yLabel = 'Percent change variance explained (%)';
        elseif strcmp(summaryMetric, 'zscoreVE')
            zscoredVE(s,:) = zscore(meanSelectedSensors(s,:));
            dataToPlot = zscoredVE;
            yl = [-3 3];
            yLabel = 'Z-scored variance explained (%)';
        end
        
        % Compute mean and standard error of variance explained across selected sensors
        figure(fH1);
        plot(range,meanSelectedSensors(s,:),'Color', lineColorSat(:,s), 'Linewidth',2); hold on;
        
        % We will also plot all individual subjects
        figure(fH2);
        subplot(2,5,s);
        
        % Plot mean with shaded error bar using 'patch' function
        se   = 100.*nanstd(thisSubjectSensorData,[],2) ./ sqrt(size(thisSubjectSensorData,2));
        ci   = 1 .* se; % use zsore=1 for 68% CI or zcore=1.95 for 95%
        lo = dataToPlot(s,:) - ci';
        hi = dataToPlot(s,:) + ci';
        
        err = patch([range, fliplr(range)], [lo, fliplr(hi)], colorCIPatch, 'FaceAlpha', 0.5, 'LineStyle',':');  hold on;
        plot(range,dataToPlot(s,:),'Color', 'r', 'Linewidth',2); hold on;
        plot([1 1], [min(yl), max(yl)], 'k');
        
        set(gca,'TickDir', 'out');
        xlabel('Scale factor');
        set(gca,'XTick', range([1 8 19]),'XTickLabel',range([1 8 19]), 'YLim', yl, 'XLim', [range(1),range(end)]);
        set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20, 'XScale', 'log'); axis square;
        title(sprintf('S%d', s));
        ylabel(yLabel);
        box off;
        
        
    end
    
    % Average/SE/CI across subjects
    averageDataToPlot = nanmean(dataToPlot,1);
    
    % Plot mean on figure with individual lines
    figure(fH1);
    plot(range,averageDataToPlot,'r','Linewidth',6);
    plot([1 1], [min(yl), max(yl)], 'k');
    
    % Add labels and make pretty
    set(gca,'TickDir', 'out');
    xlabel('Scale factor of original pRF size');
    set(gca,'XTick', range,'XTickLabel',range, 'YLim', yl, 'XLim', [range(1),range(end)]);
    set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20, 'XScale', 'log'); axis square;
    title('Variance explained by modelfit: Vary Size');
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
    
    fH3 = figure(3); clf; set(gcf,'Position',[1, 592, 838, 746]);
    plot(range,zeros(size(aveBoot)),'k','Linewidth',1); hold on;
    patch([range, fliplr(range)], [lo, fliplr(hi)],colorCIPatch, 'FaceAlpha', 0.5, 'LineStyle',':');
    plot(range,aveBoot,'r','Linewidth',5);
    plot([1 1], [min(lo), max(hi)], 'k');
    
    % Add labels and make pretty
    set(gca,'TickDir', 'out');
    xlabel('Scale factor of original pRF size');
    set(gca,'XTick', range,'XTickLabel',range, 'YLim', [min(lo), max(hi)], 'XLim', [range(1),range(end)]);
    set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20, 'XScale', 'log'); axis square;
    title('Variance explained by modelfit: Vary Size');
    ylabel(yLabel);
    box off;
    
    % Save fig
    if opt.saveFig
        
        fprintf('\n(%s): Saving figure 3 in %s\n',mfilename, saveDir);
        print(fH1, fullfile(saveDir, sprintf('fig3a_AVERAGE_varySizeSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-dpdf');
        print(fH1, fullfile(saveDir, sprintf('fig3a_AVERAGE_varySizeSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-dpng');
        
        
        print(fH2, fullfile(saveDir, sprintf('fig3a_IndividualSubjects_varySizeSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-depsc');
        print(fH2, fullfile(saveDir, sprintf('fig3a_IndividualSubjects_varySizeSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-dpng');
        
        print(fH3, fullfile(saveDir, sprintf('fig3a_BootStappedAVERAGE_varyPositionSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-dpdf');
        print(fH3, fullfile(saveDir, sprintf('fig3a_BootStappedAVERAGE_varyPositionSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-dpng');
        
        
    end
    
    return