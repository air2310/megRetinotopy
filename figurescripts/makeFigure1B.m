function makeFigure1B(dirPth,opt,saveSubDir)
% Function to create figure 1B (Predicted responses for every MEG channel)

saveDir = fullfile(dirPth.finalFig.savePth,saveSubDir,'figure1B');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% Load measured and predicted MEG responses
load(fullfile(dirPth.model.saveDataPth,opt.subfolder,'pred_resp','meanPhRefAmp10Hz'), 'meanPhRefAmp10Hz');
load(fullfile(dirPth.model.saveDataPth,opt.subfolder,'pred_resp','predMEGResponseScaled'), 'predMEGResponseScaled');
load(fullfile(dirPth.model.saveDataPth,opt.subfolder,'pred_resp','meanVarExpl.mat'),'meanVarExpl');
load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'megStimConditions.mat'), 'triggers');

% Sort the channels in the descending order of variance explained
[val, idx] = sort(meanVarExpl, 'descend');

% Select channels that have high variance explained values
nSensors = 10;

tmp = idx(~isnan(val));
topSensor=tmp(1:nSensors);
ve = val(~isnan(val));
ve_toPlot = round(ve.*100);

% axis properties
lW_orig = 3;
lW_pred = 6;
markerColor     = [0 0 0]; %[0.3010, 0.7450, 0.9330];
blinkblankColor = [0.5 0.5 0.5];

xLbl = 'Time (s)';
yLbl = 'Phase-referenced 10 Hz amplitudes (T)';
fontSize = 30;

% Define time scale
[nEpochs, nSensors] = size(meanPhRefAmp10Hz);
epochLengthSeconds = diff(opt.meg.epochStartEnd);
t = (0:nEpochs-1) .* epochLengthSeconds;

% Define blink and blank blocks
blinkIdx = (triggers.stimConditions(1:length(t))==20);
blankIdx = (triggers.stimConditions(1:length(t))==10);
blink_t = t(blinkIdx);
blank_t = t(blankIdx);

close all;

% Set predictions during blink epochs to NaN
predMEGResponseScaled(blinkIdx,:) = NaN;

% plot time series
for tt = 1:length(topSensor)
    
    % Define figure name
    figName = strcat('MEG time series: Measured vs Pred',sprintf('Sensor %d, var expl: %1.2f',topSensor(tt), ve(tt)));
    
    % Compute y limits
    tmp_yl = max(abs([min(meanPhRefAmp10Hz(:,topSensor(tt))), max(meanPhRefAmp10Hz(:,topSensor(tt)))])).*10^14;
    if (tmp_yl > 3)
        yl = [-1*tmp_yl, tmp_yl].*10^-14;
    else
        yl = [-3,3].*10^-14;
    end
    
    % Plot it!
    fH1 = figure; clf; hold all; set(gcf, 'Color', 'w', 'Position', [66,1,1855,1001], 'Name', figName); 
    
    % Plot the blank and blink periods
    for bt = 1:length(blink_t)
        patch([blink_t(bt),blink_t(bt)+epochLengthSeconds blink_t(bt)+epochLengthSeconds blink_t(bt)],[yl(1),yl(1),yl(2),yl(2)],blinkblankColor,'FaceAlpha', 0.2, 'LineStyle','none');
    end
    
    for bt = 1:length(blank_t)
        patch([blank_t(bt),blank_t(bt)+epochLengthSeconds blank_t(bt)+epochLengthSeconds blank_t(bt)],[yl(1),yl(1),yl(2),yl(2)],blinkblankColor,'FaceAlpha', 0.7, 'LineStyle','none');
    end
    
    plot(t, zeros(size(t)), 'k');
    plot(t, meanPhRefAmp10Hz(:,topSensor(tt)), 'o--','color',markerColor, 'MarkerSize',10,'MarkerEdge',markerColor,'MarkerFace',markerColor, 'LineWidth',lW_orig);
    plot(t, predMEGResponseScaled(:,topSensor(tt)), 'r', 'LineWidth',lW_pred);
    
    % Set labels, limits, legends
    xlabel(xLbl); ylabel(yLbl);
    ylim(yl); xlim([0, max(t)])
    set(gca, 'FontSize', fontSize, 'TickDir','out','TickLength',[0.010 0.010],'LineWidth',3); box off
    
    l = findobj(gca, 'Type','Line');
    legend(l([2,1]), {'Observed', 'Predicted'}, 'Location', 'NorthEastOutside', 'FontSize', 25); legend boxoff;
    
    if opt.saveFig    
        set(fH1,'Units','Inches');
        pos = get(fH1,'Position');
        set(fH1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        figurewrite(fullfile(saveDir, sprintf('MEG_time_series_Orig_Pred_sensor_%d_%d_%s',topSensor(tt),ve_toPlot(tt), opt.fNamePostFix)),[],0,'.',1);
        figurewrite(fullfile(saveDir, sprintf('MEG_time_series_Orig_Pred_sensor_%d_%d_%s',topSensor(tt),ve_toPlot(tt), opt.fNamePostFix)),[],[1 300],'.',1);
        makeFigure1B_i(topSensor(tt),saveDir, opt); % sensor location
        
    end
end

fprintf('\n(%s): Saving figure 1B in %s\n',mfilename, saveDir);

end


function makeFigure1B_i(topSensor,saveDir, opt)
% creates and saves plot showing the position of the meg sensor
% axes('Position',[0.12 0.7 0.2 0.2]);
% box on; axis off; axis image;
fH1_1 = mprfPlotHeadLayout(topSensor, false, [], false);
saveas(fH1_1, fullfile(saveDir, sprintf('sensor_location_%d%s',topSensor,opt.fNamePostFix)), 'eps');
end