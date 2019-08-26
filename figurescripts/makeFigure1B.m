function makeFigure1B(dirPth,opt)
% Function to create figure 1B (Predicted responses for every MEG channel)

varExpFile = dir(fullfile(dirPth.model.saveDataPth,'original','coherent','pred_resp','meanVarExpl.mat'));
predRespFile = dir(fullfile(dirPth.model.saveDataPth,'original','coherent','pred_resp','meanPredResponse.mat'));
origMEGData = dir(fullfile(dirPth.model.saveDataPth,'original','coherent','pred_resp','phaseReferencedMEGData.mat'));

saveSubDir = 'figure1B';
saveDir = fullfile(dirPth.finalFig.savePth,'figure1',saveSubDir);
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% check if the modelPredictions are saved in the folder. Else run
% mprf_main.m
if ~isempty(varExpFile) && ~isempty(predRespFile) && ~isempty(origMEGData)
    
    % Load measured and predicted MEG responses
    load(fullfile(origMEGData.folder,origMEGData.name),'phRefAmp10Hz');
    load(fullfile(predRespFile.folder,predRespFile.name),'meanPredResponse');
    load(fullfile(varExpFile.folder,varExpFile.name),'meanVarExpl'); 
    
    % Sort the channels in the descending order of variance explained 
    [val, idx] = sort(meanVarExpl, 'descend');
    
%%    
    % Select channels that have high variance explained values
    nSensors = 10;
%%    
    tmp = idx(~isnan(val)); topSensor=tmp(1:nSensors);
    ve = val(~isnan(val));
     
    % define time scale
    [nEpochs, ~, nSensors, ~] = size(phRefAmp10Hz);
    t = (0:nEpochs-1) .* diff(opt.meg.epochStartEnd);

    close all;
    % Calculate mean measured MEG time series from 19 runs
    meanPhRefAmp10Hz = squeeze(nanmean(phRefAmp10Hz,2));
    for tt = 1:length(topSensor)

        % Define figure properties
        % figure 
        figPos = [10 10 1920/2 1080/2]; 
        figName = strcat('MEG time series: Orig vs Pred',sprintf('Sensor %d, var expl: %1.2f',topSensor(tt), ve(tt)));
        
        % plot
        % time series
        lW_orig = 2;        
        lW_pred = 4;        
        markerColor = [0.3010, 0.7450, 0.9330];
        
        % axis properties
        %ttl = sprintf('Sensor %d, var expl: %1.2f',topSensor(tt), ve(tt));
        xLbl = 'Time (s)';
        yLbl = 'MEG response (Tesla)';
        fontSize = 20;
               
        % blink and blank blocks
        color = [0.5 0.5 0.5];
        nanIdx = find(isnan(meanPredResponse(:,topSensor(tt))));
        blinkIdx = nanIdx;
        blinkIdx(3:2:end) = blinkIdx(3:2:end) - 1;
        blinkIdx(2:2:end) = blinkIdx(2:2:end) + 1;
        blankIdx = nanIdx;
        blankIdx(1) = blinkIdx(1) + 2;
        blankIdx(3:2:end) = blinkIdx(3:2:end) + 3;
        blankIdx(2:2:end) = blinkIdx(2:2:end) + 2;
        blink_t = t(blinkIdx);
        blank_t = t(blankIdx);
               
        % Plot the figure
        fH1 = figure; set(gcf, 'Color', 'w', 'Position', figPos, 'Name', figName);hold on;
        plot(t, meanPhRefAmp10Hz(:,topSensor(tt)), 'o--','color',[0.3010, 0.7450, 0.9330], 'MarkerSize',7,'MarkerEdge',markerColor,'MarkerFace',markerColor, 'LineWidth',lW_orig);
        hold on;
        plot(t, meanPredResponse(:,topSensor(tt)), 'color',[1 0.4 0.4], 'LineWidth',lW_pred);
                
        %title(ttl);
        xlabel(xLbl); ylabel(yLbl);        
        set(gca, 'FontSize', fontSize, 'TickDir','out','LineWidth',3); box off
             
        % set x and y axis limits
        tmp_yl = max(abs([min(meanPhRefAmp10Hz(:,topSensor(tt))), max(meanPhRefAmp10Hz(:,topSensor(tt)))])).*10^14;
        if (tmp_yl > 6)
            yl = [-1*tmp_yl, tmp_yl].*10^-14;
        else
            yl = [-6,6].*10^-14;
        end
        ylim(yl); xlim([0, max(t)])
        
        hold on;
        for tmpIdx = 1:2:length(blink_t)
            patch([blink_t(tmpIdx),blink_t(tmpIdx+1) blink_t(tmpIdx+1) blink_t(tmpIdx)],[yl(1),yl(1),yl(2),yl(2)],color,'FaceAlpha', 0.2, 'LineStyle','none');
            patch([blank_t(tmpIdx),blank_t(tmpIdx+1) blank_t(tmpIdx+1) blank_t(tmpIdx)],[yl(1),yl(1),yl(2),yl(2)],color,'FaceAlpha', 0.7, 'LineStyle','none');
        end

        legend('Data', 'Prediction', 'Location', 'NorthEast'); legend boxoff;
        
        if opt.saveFig
            print(fH1, fullfile(saveDir, sprintf('MEG_time_series_Orig_Pred_sensor_%d',topSensor(tt))), '-dpng');
            makeFigure1B_i(topSensor(tt),saveDir);
        end
    end
 
    
end
 
fprintf('\n(%s): Saving figure 1B in %s\n',mfilename, saveDir);

end


function makeFigure1B_i(topSensor,saveDir)
% creates and saves plot showing the position of the meg sensor
% axes('Position',[0.12 0.7 0.2 0.2]);
% box on; axis off; axis image;
fH1_1 = mprfPlotHeadLayout(topSensor, false, [], false);
saveas(fH1_1, fullfile(saveDir, sprintf('sensor_location_%d',topSensor)), 'eps');
end