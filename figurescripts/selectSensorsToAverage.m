function sensorLoc = selectSensorsToAverage(opt, dirPth, saveDir, meanVarExpl,type)


switch type
    
    % Get sensor locations in the back
    case 'allPosterior'
        
        load(which('meg160_example_hdr.mat'))
        layout = ft_prepare_layout([],hdr);
        xpos   = layout.pos(1:157,1);
        ypos   = layout.pos(1:157,2);
        posteriorSensorLoc = (ypos<0 & xpos<1);
        sensorLoc = find(posteriorSensorLoc);
        
    % Get sensor locations in the back with positive values
    case 'allPosteriorPositive'
        % get sensor locations
        load(which('meg160_example_hdr.mat'))
        layout = ft_prepare_layout([],hdr);
        xpos   = layout.pos(1:157,1);
        ypos   = layout.pos(1:157,2);
        
        posteriorSensorLoc = (ypos<0 & xpos<1);
        positiveSensors    = (meanVarExpl>0);
        meanVarExpl(~positiveSensors) = NaN;
        sensorLoc  = union(positiveSensors, repmat(posteriorSensorLoc', [size(meanVarExpl,1), 1]), 'rows');
        sensorLoc = find(sensorLoc==1);
    % Selecting union of top 10 sensors from all iterations with positive values
    case 'top10Positive'
        
        positiveSensors = (meanVarExpl>0);
        meanVarExpl(~positiveSensors)=NaN;
       [val,idx] = sort(meanVarExpl,2,'descend');
       nanMask = isnan(val);
       for ii = 1:size(val,1)
           countStart = sum(nanMask(ii,:));
           if countStart > (size(meanVarExpl,2)-10)
               sensorLoc(ii,:) = [idx(ii,countStart+1:size(meanVarExpl,2)), NaN(1,10-diff([countStart,size(meanVarExpl,2)]))];
           else
               sensorLoc(ii,:) = idx(ii,countStart+1:countStart+10);
           end
       end
    % Selecting union of top 10 sensors from all iterations
    case 'top10'
        meanVarExpl(isnan(meanVarExpl))=0;
        [~,idx] = sort(meanVarExpl,2,'descend');
        sensorLoc = unique(idx(:,1:10));
        
    % Selecting union of top 10 sensors from all iterations
    case 'top5'
        meanVarExpl(isnan(meanVarExpl))=0;
        [~,idx] = sort(meanVarExpl,2,'descend');
        sensorLoc = unique(idx(:,1:5));
        
    % Get top 10 sensors from splithalf reliability
    case 'reliable10'
        load(fullfile(dirPth.model.saveDataPth, 'splitHalfAmpReliability1000'), 'splitHalfAmpCorrelation');
        splitHalfAmpCorrelation(isnan(splitHalfAmpCorrelation))=0;
        
        [~,idx] = sort(splitHalfAmpCorrelation, 'descend');
        sensorLoc = unique(idx(1:10));
        
end

if size(sensorLoc,2)>1
    for ii = 1:size(sensorLoc,1)
        fh = mprfPlotHeadLayout(sensorLoc(ii,~isnan(sensorLoc(ii,:))),0,[]);
        
        if opt.saveFig
            print(gcf, fullfile(saveDir, sprintf('fig2c_%s_varyPositionSensors%s_%s_%d', dirPth.subjID, opt.fNamePostFix, type, ii)), '-dpdf');
        end
        close(fh)
    end
else 
    fh = mprfPlotHeadLayout(sensorLoc',0,[]);
    if opt.saveFig
        print(gcf, fullfile(saveDir, sprintf('fig2c_%s_varyPositionSensors%s_%s', dirPth.subjID, opt.fNamePostFix, type)), '-dpdf');
    end
    
    close(fh)
end

return
