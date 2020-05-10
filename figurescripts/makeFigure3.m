function makeFigure3(dirPth, opt, sensorsToAverage)
% Function to make that plots the effect of scaling the original pRF
% sizes (sigma).

if opt.saveFig
    [pth, ~] = fileparts(dirPth.model.saveFigPth);
    saveSubDir = ['figure3_' opt.subfolder];
    saveDir = fullfile(pth, 'finalfig', saveSubDir);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
end

% Define plotting params
color = [0.5 0.5 0.5];
yl    = [-10 50];
clim  = [0 50];
interpmethod = 'nearest'; % can also be [] for 'v4' --> smooth interpolation

% Load variance explained
load(fullfile(dirPth.model.saveDataPth, opt.subfolder, 'pred_resp', 'meanVarExpl'), 'meanVarExpl');

% Get range of permutations
range   = opt.vary.size;

% What sensors are we averaging?
sensorLoc = selectSensorsToAverage(opt, dirPth, saveDir, meanVarExpl, sensorsToAverage);

% Plot data for sensors over the back of the head
if strcmp(sensorsToAverage, 'top10Positive')
    for ii = 1:size(sensorLoc,1)
        curSensorLoc = sensorLoc(ii,:);
        dataToPlot(ii,1:sum(~isnan(curSensorLoc))) = meanVarExpl(ii,curSensorLoc(~isnan(curSensorLoc)));
    end
else
    dataToPlot = meanVarExpl(:,sensorLoc);
end

% Compute mean and standard error of variance explained across selected sensors
mean_varexpl = nanmean(dataToPlot,2);
se_varexpl   = nanstd(dataToPlot,0,2) ./ sqrt(size(dataToPlot,2));
ci_varexpl   = 1 .* se_varexpl;  % use zsore=1 for 68% CI or zcore=1.95 for 95%
lo           = 100.*(mean_varexpl - ci_varexpl);
hi           = 100.*(mean_varexpl + ci_varexpl);

% Plot mean with shaded error bar using 'patch' function
fH1 = figure(1); clf; set(fH1, 'Color', 'w', 'Position', [66,1,1855,1001], 'Name', 'Vary pRF size');
err = patch([range, fliplr(range)], [lo', fliplr(hi')], color, 'FaceAlpha', 0.5, 'LineStyle',':');
hold on;
plot(range,100.*mean_varexpl,'r','Linewidth',3);
plot(range, zeros(size(range)), 'k')
plot([1 1], yl, 'k')

% Add labels and make pretty
set(gca,'TickDir', 'out');
xlabel('Scale factor');
set(gca,'XTick', range,'XTickLabel',range, 'YLim', yl, 'XLim', [range(1),range(end)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20, 'XScale', 'log'); axis square;
title('Variance explained by modelfit: Vary Size');
ylabel('Variance explained (%)');


%% Plot meshes
rows = 3;
cols = round(length(range)/rows)+1;
fH2 = figure(2); clf; set(fH2, 'Color', 'w', 'Position', [ 136, 96, 2000,  1138],  'Name', 'Vary pRF size'); hold all;

for ii = 1:length(range)
    % Select data
    meshDataToPlot = meanVarExpl(ii,:).*100;
    
    % Get subplot
    subplot(rows, cols, ii);
    
    megPlotMap(meshDataToPlot,clim,fH2,'parula',...
        sprintf('%1.2fx',range(ii)),[],[],'interpmethod',interpmethod);
    c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
    pos = c.Position; set(c, 'Position', [pos(1)+0.03 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])   
end

%% Save figures if requestsed
if opt.saveFig
    fprintf('\n(%s): Saving figure 3 in %s\n',mfilename, saveDir);
    print(fH1, fullfile(saveDir, sprintf('fig3a_%s_varySizeSummary%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), '-dpdf');
    print(fH2, fullfile(saveDir, sprintf('fig3b_%s_varySizeMeshes%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), '-dpng');
end


return