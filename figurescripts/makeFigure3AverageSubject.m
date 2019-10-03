function makeFigure3AverageSubject(opt, sensorsToAverage)
% Function to make average across subjects for Figure 2 from manuscript, 
% plotting variance explained by the model as a function of polar angle 
% rotations around the fovea of the original estimated pRF centers.

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
fH1   = figure(1); clf; set(fH1, 'Color', 'w', 'Position', [1000, 592, 838, 746]);
colorLine = [0.8 0.8 0.8];
colorPatch = [0.5 0.5 0.5];

for s = 1:length(subjects)
    
    subjectID = subjects{s};
    
    dirPth = loadPaths(subjectID);
    
    % Load variance explained file
    load(fullfile(dirPth.model.saveDataPth, 'vary_size','coherent','pred_resp', 'meanVarExpl'));
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
    dataToPlot = squeeze(varexpl(s,:,sensorLoc{s}));
    
    % Compute mean and standard error of variance explained across selected sensors
    meanSelectedSensors(s,:) = nanmean(dataToPlot,2);
    plot(range,meanSelectedSensors(s,:),'Color', colorLine, 'Linewidth',2); hold on;
    
end

% Average/SE/CI across subjects
averageSubjectVarExpl = nanmean(meanSelectedSensors,1);
averageSubjectSE = nanstd(meanSelectedSensors) ./ sqrt(length(subjects));
averageSubjectCI = 1.96.* averageSubjectSE;


% % Plot shaded error bar using 'patch' function
% lo = averageSubjectVarExpl - averageSubjectCI;
% hi = averageSubjectVarExpl + averageSubjectCI;
% 
% err = patch([range, fliplr(range)], [lo', fliplr(hi')], colorPatch, 'FaceAlpha', 0.5, 'LineStyle',':');

% Plot mean
plot(range,averageSubjectVarExpl,'r','Linewidth',4);

% Add labels and make pretty
set(gca,'TickDir', 'out');
xlabel('Position (deg)');
set(gca,'XTick', range,'XTickLabel',range, 'YLim', [0 inf], 'XLim', [range(1),range(end)]);
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20, 'XScale', 'log'); axis square;
title('Variance explained by modelfit: Vary Size');
ylabel('Variance explained (%)');
box off;

% Save fig
if opt.saveFig
    saveDir = fullfile(dirPth.finalFig.savePthAverage, 'figure3');
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    fprintf('\n(%s): Saving figure 3 in %s\n',mfilename, saveDir);

    print(fH1, fullfile(saveDir, sprintf('fig3a_AVERAGE_varyPositionSummary%s_%s', opt.fNamePostFix, sensorsToAverage)), '-dpng');
end

return