function makeFigure2AverageSubject(opt, sensorsToAverage, summaryMetric)
% Function to make average across subjects for Figure 2 from manuscript, 
% plotting variance explained by the model as a function of polar angle 
% rotations around the fovea of the original estimated pRF centers.

% check inputs
if ~exist('summaryMetric', 'var') || isempty(summaryMetric)
    summaryMetric = 'meanVE';
end

subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
            'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

% Define the range of rotations
range   = opt.vary.position;
        
% Allocate space
varexpl = NaN(length(subjects),length(range), 157);
sensorLoc = cell(length(subjects),1);

% Get sensor locations in the back
load(which('meg160_example_hdr.mat'))
layout = ft_prepare_layout([],hdr);
xpos = layout.pos(1:157,1);
ypos = layout.pos(1:157,2);

% Figure for all subjects in one plot
fH1   = figure(1); clf; set(fH1, 'Color', 'w', 'Position', [1000, 592, 838, 746]);

% Figure for all subjects plotted separately
if  strcmp(summaryMetric, 'meanVE')  
    sz = get(0, 'screensize');
    fH2 = figure; set(gcf,'Position', sz);
    colorCIPatch = [0.5 0.5 0.5];
end

lineColorSat = repmat(linspace(0.3,0.9,10), [3 1]);

for s = 1:length(subjects)
    
    subjectID = subjects{s};
    
    dirPth = loadPaths(subjectID);
    
    % Load variance explained file
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp', 'meanVarExpl'));
    varexpl(s,:,:) = meanVarExpl;
    
    if strcmp(sensorsToAverage, 'top10')
        % Get top 10 sensors
        tmp = squeeze(varexpl(s,:,:));
        tmp(isnan(tmp))=0;
        [~,idx] = sort(tmp,2,'descend');
        sensorLoc{s} = unique(idx(:,1:10)); % selecting union of top 10 sensors from all iterations
    elseif strcmp(sensorsToAverage, 'allPosterior')
        sensorLoc{s} = (ypos<0 & xpos<1);
    end
    
    % Plot data for sensors over the back of the head
    thisSubjectSensorData = squeeze(varexpl(s,:,sensorLoc{s}));
    
    % Compute summary metrics of variance explained across selected sensors
    meanSelectedSensors(s,:) = 100*nanmean(thisSubjectSensorData,2);

    if strcmp(summaryMetric, 'meanVE')
        dataToPlot = meanSelectedSensors;      
        yl = [0 45];
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
    
    figure(fH1);
    plot(range,dataToPlot(s,:),'Color', lineColorSat(:,s), 'Linewidth',2); hold on;
    
    % If we plot mean VE, then also make plot with all individual subjects
    if  strcmp(summaryMetric, 'meanVE')    
        figure(fH2);
        subplot(2,5,s);

        % Plot mean with shaded error bar using 'patch' function
        se   = 100.*nanstd(thisSubjectSensorData,0,2) ./ sqrt(size(thisSubjectSensorData,2));
        ci   = 1.96 .* se;      
        lo = dataToPlot(s,:) - ci';
        hi = dataToPlot(s,:) + ci';

        err = patch([range, fliplr(range)], [lo, fliplr(hi)], colorCIPatch, 'FaceAlpha', 0.5, 'LineStyle',':');  hold on;
        plot(range,dataToPlot(s,:),'Color', 'r', 'Linewidth',2); hold on;

        set(gca,'TickDir', 'out'); xlabel('Position (deg)');
        set(gca,'XTick', range([1 5 9]),'XTickLabel',rad2deg(range([1 5 9])), 'YLim', yl, 'XLim', [range(1),range(end)]);
        set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
        title(sprintf('S%d', s));
        ylabel(yLabel); box off;
    end
    
end

% Average/SE/CI across subjects
averageDataToPlot          = nanmean(dataToPlot,1);

% averageSubjectSE = nanstd(meanSelectedSensors) ./ sqrt(length(subjects));
% averageSubjectCI = 1.96.* averageSubjectSE;

% % Plot shaded error bar using 'patch' function
% lo = averageSubjectVarExpl - averageSubjectCI;
% hi = averageSubjectVarExpl + averageSubjectCI;
% 
% err = patch([range, fliplr(range)], [lo', fliplr(hi')], colorPatch, 'FaceAlpha', 0.5, 'LineStyle',':');

% Plot mean
figure(fH1);
plot(range,averageDataToPlot,'r','Linewidth',5);

% Add labels and make pretty
set(gca,'TickDir', 'out');
xlabel('Position (deg)');
set(gca,'XTick', range,'XTickLabel',rad2deg(range), 'YLim', yl, 'XLim', [range(1),range(end)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
title('Variance explained by modelfit: Vary Position');
ylabel(yLabel);
box off;

% Save fig
if opt.saveFig
    saveDir = fullfile(dirPth.finalFig.savePthAverage, 'figure2');
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    fprintf('\n(%s): Saving figure 2 in %s\n',mfilename, saveDir);

    print(fH1, fullfile(saveDir, sprintf('fig2a_AVERAGE_varyPositionSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-depsc');
    print(fH1, fullfile(saveDir, sprintf('fig2a_AVERAGE_varyPositionSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-dpng');

    if  strcmp(summaryMetric, 'meanVE')    
        print(fH2, fullfile(saveDir, sprintf('fig2a_IndividualSubjects_varyPositionSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-depsc');
        print(fH2, fullfile(saveDir, sprintf('fig2a_IndividualSubjects_varyPositionSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-dpng');
    end
end



return