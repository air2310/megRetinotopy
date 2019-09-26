function makeFigure3(dirPth, opt,sensorsToAverage)
% Function to make that plots the effect of scaling the original pRF
% sizes (sigma).

varexpl = load(fullfile(dirPth.model.saveDataPth, 'vary_size','coherent','pred_resp', 'meanVarExpl'));
varexpl = varexpl.meanVarExpl;

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
    tmp = varexpl;
    tmp(isnan(tmp))=0;
    [val,idx] = sort(tmp,2,'descend');
    sensorLoc = unique(idx(:,1:10)); % selecting union of top 10 sensors from all iterations
end

% Plot data for sensors over the back of the head
dataToPlot = varexpl(:,sensorLoc);

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
set(gca,'XTick', range,'XTickLabel',range, 'YLim', [0 0.25], 'XLim', [range(1),range(end)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20, 'XScale', 'log'); axis square;
title('Variance explained by modelfit: Vary Size');
ylabel('Variance explained (%)');



%% Plot meshes
rows = 3;
cols = round(length(range)/rows)+1;
fH2 = figure(2); clf; set(fH2, 'Color', 'w', 'Position', [ 136, 96, 2000,  1138]); hold all;

clim = max(varexpl(range==1,:));
interpmethod = 'nearest'; % can also be 'v4' for smooth interpolation

for ii = 1:length(range)
    % Select data
    meshDataToPlot = varexpl(ii,:);
    
    % Get subplot
    subplot(rows, cols, ii);
    
    megPlotMap(meshDataToPlot,[0 0.6],fH2,'parula',...
        sprintf('%1.2fx',range(ii)),[],[],'interpmethod',interpmethod);
end


%% Plot the sensors selected for averaging
fH3 = mprfPlotHeadLayout(sensorLoc',0,[]);

%%
if opt.saveFig
    [pth, folder] = fileparts(dirPth.model.saveFigPth);
    saveDir = fullfile(pth, 'finalfig', 'figure3');
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    fprintf('\n(%s): Saving figure 3 in %s\n',mfilename, saveDir);
    
    print(fH1, fullfile(saveDir, sprintf('fig3a_%s_varySizeSummary%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), '-dpng');
    print(fH2, fullfile(saveDir, sprintf('fig3b_%s_varySizeMeshes%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), '-dpng');
    print(fH3, fullfile(saveDir, sprintf('fig3b_%s_varySizeSensors%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)), '-dpng');


end


return