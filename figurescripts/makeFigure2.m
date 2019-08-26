function makeFigure2(dirPth, opt)
% Function to make that plots the effect of rotating the original pRF
% polar angle positions around the fovea.

varexpl = load(fullfile(dirPth.model.saveDataPth, 'vary_position','coherent','pred_resp', 'meanVarExpl'));
varexpl = varexpl.meanVarExpl;

range   = opt.vary.position;

%% Plot summary
% Get sensor locations in the back
load(which('meg160_example_hdr.mat'))
layout = ft_prepare_layout([],hdr);
xpos = layout.pos(1:157,1);
ypos = layout.pos(1:157,2);
sensorLocBack = (ypos<0 & xpos<1);

% Plot data for sensors over the back of the head
dataToPlot = varexpl(:,sensorLocBack);

% Compute mean and standard error of variance explained across selected sensors
mean_varexpl = nanmean(dataToPlot,2);
se_varexpl   = nanstd(dataToPlot,0,2) ./ sqrt(size(dataToPlot,2));
ci_varexpl   = 1.96 .* se_varexpl;

fH1 = figure(1); clf; set(fH1, 'Color', 'w', 'Position', [1000, 592, 838, 746]);

% Plot mean with shaded error bar using 'patch' function
lo = mean_varexpl - ci_varexpl;
hi = mean_varexpl + ci_varexpl;
color = [0.5 0.5 0.5];
err = patch([range, fliplr(range)], [lo', fliplr(hi')], color, 'FaceAlpha', 0.5, 'LineStyle',':');
hold on;
plot(range,mean_varexpl,'r','Linewidth',3);

% Add labels and make pretty
set(gca,'TickDir', 'out');
xlabel('Position (deg)');
set(gca,'XTick', range,'XTickLabel',rad2deg(range), 'YLim', [0 0.25], 'XLim', [range(1),range(end)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
title('Variance explained by modelfit: Vary Position');
ylabel('Variance explained (%)');



%% Plot meshes
rows = 2;
cols = round(length(range)/rows);
fH2 = figure(2); clf; set(fH2, 'Color', 'w', 'Position', [326,584,1234,754]); hold all;

clim = max(varexpl(range==0,:));
interpmethod = 'nearest'; % can also be 'v4' for smooth interpolation

for ii = 1:length(range)
    % Select data
    meshDataToPlot = varexpl(ii,:);
    
    % Get subplot
    subplot(rows, cols, ii);
    
    megPlotMap(meshDataToPlot,[0 clim],fH2,'parula',...
        range(ii),[],[],'interpmethod',interpmethod);
end


if opt.saveFig
    if ~exist(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'msFigs'), 'dir')
        mkdir(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'msFigs')); end
    
    print(fH1, fullfile(dirPth.model.saveFigPth, opt.subfolder, 'msFigs', sprintf('fig2a_%s_varyPositionSummary%s', dirPth.subjID, opt.fNamePostFix)), '-dpng');
    print(fH2, fullfile(dirPth.model.saveFigPth, opt.subfolder, 'msFigs', sprintf('fig2b_%s_varyPositionMeshes%s', dirPth.subjID, opt.fNamePostFix)), '-dpng');

end


return