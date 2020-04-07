function makeFigure2(dirPth, opt, sensorsToAverage)
% Function to make Figure 2 from manuscript, plotting variance explained by
% the model as a function of polar angle rotations around the fovea of the 
% original estimated pRF centers.

if opt.saveFig
    [pth, ~] = fileparts(dirPth.model.saveFigPth);
    saveSubDir = ['figure2_' opt.subfolder];
    saveDir = fullfile(pth, 'finalfig', saveSubDir);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
end

% Load variance explained file
load(fullfile(dirPth.model.saveDataPth, opt.subfolder, 'pred_resp', 'meanVarExpl'), 'meanVarExpl');

% Define the range of rotations
range   = opt.vary.position;

% What sensors are we averaging?
if strcmp(sensorsToAverage, 'allPosterior')
    % Get sensor locations in the back
    load(which('meg160_example_hdr.mat'))
    layout = ft_prepare_layout([],hdr);
    xpos = layout.pos(1:157,1);
    ypos = layout.pos(1:157,2);
    sensorLoc = (ypos<0 & xpos<1);
elseif strcmp(sensorsToAverage, 'top10')
    % Get top 10 sensors
    tmp = meanVarExpl;
    tmp(isnan(tmp))=0;
    [~,idx] = sort(tmp,2,'descend');
    sensorLoc = unique(idx(:,1:10)); % selecting union of top 10 sensors from all iterations
end

% Plot data for sensors over the back of the head
dataToPlot = meanVarExpl(:,sensorLoc);

% Compute mean and standard error of variance explained across selected sensors
mean_varexpl = nanmean(dataToPlot,2);
se_varexpl   = nanstd(dataToPlot,0,2) ./ sqrt(size(dataToPlot,2));
ci_varexpl   = 1.96 .* se_varexpl;

% Plot mean with shaded error bar using 'patch' function
lo = 100.*(mean_varexpl - ci_varexpl);
hi = 100.*(mean_varexpl + ci_varexpl);
color = [0.5 0.5 0.5];

fH1 = figure(1); clf; set(fH1, 'Color', 'w', 'Position', [66,1,1855,1001], 'Name', 'Vary pRF position');
err = patch([range, fliplr(range)], [lo', fliplr(hi')], color, 'FaceAlpha', 0.5, 'LineStyle',':');
hold all;
plot(range,100.*mean_varexpl,'r','Linewidth',3);
plot(range, zeros(size(range)), 'k')
% Add labels and make pretty
yl = [0 50];
if max(100.*mean_varexpl)>yl(2)
    yl = [0 max(100.*mean_varexpl)+5];
end
if min(100.*mean_varexpl)<yl(1)
    yl = [min(100.*mean_varexpl)-5 yl(2)];
end
set(gca,'TickDir', 'out');
xlabel('Position (deg)');
set(gca,'XTick', range,'XTickLabel',rad2deg(range), 'YLim', yl, 'XLim', [range(1),range(end)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
title('Variance explained by modelfit: Vary Position');
ylabel('Variance explained (%)');


%% Plot meshes
rows = 2;
cols = round(length(range)/rows);
fH2 = figure(2); clf; set(fH2, 'Color', 'w', 'Position', [326,584,1234,754], 'Name', 'Vary pRF position'); hold all;

clim = yl;
%interpmethod = 'nearest'; % can also be [] for 'v4' --> smooth interpolation
interpmethod = []; % using the default 'v4' interpolation

for ii = 1:length(range)
    % Select data
    meshDataToPlot = meanVarExpl(ii,:).*100;
        
    % Get subplot
    fH2 = figure(2); hold all;
    subplot(rows, cols, ii);
    
    megPlotMap(meshDataToPlot,clim,fH2,'parula',...
        range(ii),[],[],'interpmethod',interpmethod);
    c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
    pos = c.Position; set(c, 'Position', [pos(1)+0.04 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])

end


%% Plot the sensors selected for averaging
if strcmp(sensorsToAverage, 'top10')
    fH3 = mprfPlotHeadLayout(sensorLoc',0,[]);
end

if opt.saveFig

    fprintf('\n(%s): Saving figure 2 in %s\n',mfilename, saveDir);
    
%     set(fH1,'Units','Inches');
%     pos = get(fH1,'Position');
%     set(fH1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(fH1, fullfile(saveDir, sprintf('fig2a_%s_varyPositionSummary%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), '-dpdf');
    
    figure(fH1);
    figurewrite(fullfile(saveDir, sprintf('fig2a_%s_varyPositionSummary%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), [],[1 300],'.',1);
    figurewrite(fullfile(saveDir, sprintf('fig2a_%s_varyPositionSummary%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), [],0,'.',1);

    figure(fH2);
    figurewrite(fullfile(saveDir, sprintf('fig2b_%s_varyPositionMeshes%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)),[],[1 300],'.',1);
    figurewrite(fullfile(saveDir, sprintf('fig2b_%s_varyPositionMeshes%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)),[],0,'.',1);
    fprintf('\n saving figure 1D in %s',saveDir);
    
    if strcmp(sensorsToAverage, 'top10')
        print(fH3, fullfile(saveDir, sprintf('fig2c_%s_varyPositionSensors%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), '-dpdf');
    end
end

return