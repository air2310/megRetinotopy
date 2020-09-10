function [] = makeSupplementFigureS5(sensorsToAverage, summaryMetric, opt)
% Function to plot Supplement Figure S5, showing average variance explained 
% by the model as a function of scale factor of the original estimated pRF 
% size, for individual subjects.

% Check inputs
if ~exist('summaryMetric', 'var') || isempty(summaryMetric)
    summaryMetric = 'meanVE';
end

if ~exist('sensorsToAverage', 'var') || isempty(sensorsToAverage)
    sensorsToAverage = 'top10';
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

% Figure for all subjects plotted separately
sz = get(0, 'screensize');
fH1 = figure; set(gcf,'Position', sz);
colorCIPatch = [0.5 0.5 0.5];


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
    
    % We will also plot all individual subjects
    figure(fH1);
    subplot(2,5,s);
    
    % Plot mean with shaded error bar using 'patch' function
    se   = 100.*nanstd(thisSubjectSensorData,[],2) ./ sqrt(size(thisSubjectSensorData,2));
    ci   = 1 .* se; % use zsore=1 for 68% CI or zcore=1.95 for 95%
    lo = dataToPlot(s,:) - ci';
    hi = dataToPlot(s,:) + ci';
    
    err = patch([range, fliplr(range)], [lo, fliplr(hi)], colorCIPatch, 'FaceAlpha', 0.5, 'LineStyle',':');  hold on;
    plot(range,dataToPlot(s,:),'Color', 'r', 'Linewidth',2); hold on;
    plot(range, zeros(size(dataToPlot(s,:))), 'k', 'LineWidth', 1);
    plot([1 1], [min(yl), max(yl)], 'k');
    
    set(gca,'TickDir', 'out');
    xlabel('Scale factor');
    set(gca,'XTick', range([1 8 19]),'XTickLabel',range([1 8 19]), 'YLim', yl, 'XLim', [range(1),range(end)]);
    set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20, 'XScale', 'log'); axis square;
    title(sprintf('S%d', s));
    ylabel(yLabel);
    box off;
  
end % subjects

% Save figures
if opt.saveFig
    % Define folder to save figures
    saveSubDir = ['SupplFigure4_' opt.subfolder];
    saveDir = fullfile(dirPth.finalFig.savePthAverage, saveSubDir);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    
    fprintf('\n(%s): Saving Supplement Figure 4 in %s\n',mfilename, saveDir);
  
    print(fH1, fullfile(saveDir, sprintf('SuppFigure4_IndividualSubjects_varySizeSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-depsc');
    print(fH1, fullfile(saveDir, sprintf('SuppFigure4_IndividualSubjects_varySizeSummary%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-dpng');
    
end

return