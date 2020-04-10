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

% Load variance explained
load(fullfile(dirPth.model.saveDataPth, opt.subfolder, 'pred_resp', 'meanVarExpl'), 'meanVarExpl');

% Get range of permutations
range   = opt.vary.size;

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
    [val,idx] = sort(tmp,2,'descend');
    sensorLoc = unique(idx(:,1:10)); % selecting union of top 10 sensors from all iterations
end

% Plot data for sensors over the back of the head
dataToPlot = meanVarExpl(:,sensorLoc);

% Compute mean and standard error of variance explained across selected sensors
mean_varexpl = nanmean(dataToPlot,2);
se_varexpl   = nanstd(dataToPlot,0,2) ./ sqrt(size(dataToPlot,2));
ci_varexpl   = 1 .* se_varexpl;  % use zsore=1 for 68% CI or zcore=1.95 for 95%


fH1 = figure(1); clf; set(fH1, 'Color', 'w', 'Position', [66,1,1855,1001], 'Name', 'Vary pRF size');

% Plot mean with shaded error bar using 'patch' function
lo = 100.*(mean_varexpl - ci_varexpl);
hi = 100.*(mean_varexpl + ci_varexpl);
color = [0.5 0.5 0.5];
err = patch([range, fliplr(range)], [lo', fliplr(hi')], color, 'FaceAlpha', 0.5, 'LineStyle',':');
hold on;
plot(range,100.*mean_varexpl,'r','Linewidth',3);

% Add labels and make pretty
yl = [0 50];
if max(100.*mean_varexpl)>yl(2)
    yl = [0 max(100.*mean_varexpl)+5];
end
if min(100.*mean_varexpl)<yl(1)
    yl = [min(100.*mean_varexpl)-5 yl(2)];
end

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

clim = [0 50];
%interpmethod = 'nearest'; % can also be 'v4' for smooth interpolation
interpmethod = []; % using the default 'v4' interpolation

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


%% Plot the sensors selected for averaging

if strcmp(sensorsToAverage, 'top10')
    fH3 = mprfPlotHeadLayout(sensorLoc',0,[]);
end

%%
if opt.saveFig

    fprintf('\n(%s): Saving figure 3 in %s\n',mfilename, saveDir);
    
%     set(fH1,'Units','Inches');
%     pos = get(fH1,'Position');
%     set(fH1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
    print(fH1, fullfile(saveDir, sprintf('fig3a_%s_varySizeSummary%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), '-dpdf');
    print(fH2, fullfile(saveDir, sprintf('fig3b_%s_varySizeMeshes%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), '-dpng');
    if strcmp(sensorsToAverage, 'top10')
        print(fH3, fullfile(saveDir, sprintf('fig3b_%s_varySizeSensors%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), '-dpdf');
    end

end


return