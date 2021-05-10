function makeFigure5(dirPth, opt, sensorsToAverage)
% Function to make Figure 2 from manuscript, plotting variance explained by
% the model as a function of polar angle rotations around the fovea of the 
% original estimated pRF centers.

if opt.saveFig
    [pth, ~] = fileparts(dirPth.model.saveFigPth);
    saveSubDir = ['Figure5_' opt.subfolder];
    saveDir = fullfile(pth, 'finalfig', saveSubDir);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
end

% Define plotting params
color = [0.5 0.5 0.5];
yl    = [-5 50];

% Load variance explained file
load(fullfile(dirPth.model.saveDataPth, opt.subfolder, 'pred_resp', 'meanVarExpl'), 'meanVarExpl');

% Define the range of rotations
range   = opt.vary.position;

% What sensors are we averaging?
sensorLoc = selectSensorsToAverage(opt, dirPth, saveDir, meanVarExpl, sensorsToAverage);

% Select data to plot
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
ci_varexpl   = 1 .* se_varexpl;  % use zsore=1 for 68% CI or zcore=1.96 for 95%
lo           = 100.*(mean_varexpl - ci_varexpl);
hi           = 100.*(mean_varexpl + ci_varexpl);

% Plot mean with shaded error bar using 'patch' function
fH1 = figure(1); clf; set(fH1, 'Color', 'w', 'Position', [66,1,1855,1001], 'Name', 'Vary pRF position');
err = patch([range, fliplr(range)], [lo', fliplr(hi')], color, 'FaceAlpha', 0.5, 'LineStyle',':');
hold all;
plot(range,100.*mean_varexpl,'r','Linewidth',3);
plot(range, zeros(size(range)), 'k')
plot([0 0 ], yl, 'k')

% Add labels and make pretty
set(gca,'TickDir', 'out');
xlabel('Rotation angle (deg)');
set(gca,'XTick', range,'XTickLabel',rad2deg(range), 'YLim', yl, 'XLim', [range(1),range(end)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
title('Variance explained by modelfit: Vary Position');
ylabel('Variance explained (%)');


%% Save figures if requested
if opt.saveFig
    [pth, ~] = fileparts(dirPth.model.saveFigPth);
    saveSubDir = ['Figure5_' opt.subfolder];
    saveDir = fullfile(pth, 'finalfig', saveSubDir, sensorsToAverage);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end

    fprintf('\n(%s): Saving Figure 5 in %s\n',mfilename, saveDir);
    print(fH1, fullfile(saveDir, sprintf('Figure5A_%s_varyPositionSummary%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), '-dpdf');
    figurewrite(fullfile(saveDir, sprintf('Figure5A_%s_varyPositionSummary%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), [],[1 300],'.',1);
%     figurewrite(fullfile(saveDir, sprintf('Figure5A_%s_varyPositionSummary%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), [],0,'.',1);

end

return