function makeFigure1B(dirPth,opt)
% Function to create figure 1B (Predicted responses for every MEG channel)

varExpFile = dir(fullfile(dirPth.model.saveDataPth,'original','pred_resp','meanVarExpl.mat'));
predRespFile = dir(fullfile(dirPth.model.saveDataPth,'original','pred_resp','meanPredResponse.mat'));
origMEGData = dir(fullfile(dirPth.model.saveDataPth,'original','pred_resp','phaseReferencedMEGData.mat'));

saveSubDir = 'figure1B';
saveDir = fullfile(dirPth.finalFig.savePth,'figure1',saveSubDir);
if ~exist(saveSubDir,'dir')
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
    t = (0:nEpochs-1) .* diff(opt.epochStartEnd);

    close all;
    % Calculate mean measured MEG time series from 19 runs
    meanPhRefAmp10Hz = squeeze(nanmean(phRefAmp10Hz,2));
    for tt = 1:length(topSensor)
        
        fH1 = figure; set(gcf, 'Color', 'w', 'Position', [10 10 1920/2 1080/2], 'Name', 'MEG time series: Orig vs Pred');

        plot(t, meanPhRefAmp10Hz(:,topSensor(tt)), 'ko-', 'LineWidth',2);
        hold on;
        plot(t, meanPredResponse(:,topSensor(tt)), 'r', 'LineWidth',4);
        
        title(sprintf('Sensor %d, var expl: %1.2f',topSensor(tt), ve(tt)));
        xlabel('Time (s)'); ylabel('MEG response (Tesla)');
        
        set(gca, 'FontSize', 14, 'TickDir','out'); box off
        
        tmp_yl = max(abs([min(meanPhRefAmp10Hz(:,topSensor(tt))), max(meanPhRefAmp10Hz(:,topSensor(tt)))])).*10^14;
        if (tmp_yl > 6)
            yl = [-1*tmp_yl, tmp_yl].*10^-14; 
        else
            yl = [-6,6].*10^-14;
        end
        ylim(yl); xlim([0, max(t)])
        
        legend('Data', 'Prediction', 'Location', 'NorthEast'); legend boxoff;

        % create a inset plot showing the position of the meg sensor
        axInset = axes('Position',[0.2 0.7 0.2 0.2]);
        box on; axis off; axis image;
        holdFig.flag = 1; holdFig.figh = fH1; holdFig.c = 'r';
        mprfPlotHeadLayout(topSensor(tt), false, [], false,holdFig);
        
        print(fH1, fullfile(saveDir, sprintf('MEG_time_series_Orig_Pred_sensor_%d',topSensor(tt))), '-dpng');
    end
 
    
end
 
end