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
       
    % Selecting union of top 15 sensors from all iterations
    case 'top15'
        meanVarExpl(isnan(meanVarExpl))=0;
        [~,idx] = sort(meanVarExpl,2,'descend');
        sensorLoc = unique(idx(:,1:15))';
       
    % Selecting union of top 10 sensors from all iterations
    case 'top10'
        meanVarExpl(isnan(meanVarExpl))=0;
        [~,idx] = sort(meanVarExpl,2,'descend');
        sensorLoc = unique(idx(:,1:10))';
        
    % Selecting union of top 10 sensors from all iterations
    case 'top5'
        meanVarExpl(isnan(meanVarExpl))=0;
        [~,idx] = sort(meanVarExpl,2,'descend');
        sensorLoc = unique(idx(:,1:5))';
        
    % Get top 10 sensors from splithalf reliability
    case 'top10reliable'
        load(fullfile(dirPth.model.saveDataPth, 'splitHalfAmpReliability1000'), 'splitHalfAmpCorrelation');
        if size(splitHalfAmpCorrelation,1)>1
            splitHalfAmpCorrelation = mean(splitHalfAmpCorrelation,1,'omitnan');
        end
        splitHalfAmpCorrelation(isnan(splitHalfAmpCorrelation))=0;
        [~,idx] = sort(splitHalfAmpCorrelation, 'descend');
        sensorLoc = unique(idx(1:10));
        
    % Get top 10 sensors from fit using original pRF params
    case 'top10orig'
        if size(meanVarExpl,1) == 9 
            origIdx = 5;
        else
            origIdx = 8;
        end
         meanVarExpl(isnan(meanVarExpl))=0;
        [~,idx] = sort(meanVarExpl(origIdx,:),'descend');
        sensorLoc = unique(idx(1:10));
end

if size(sensorLoc,1)>1
    for ii = 1:size(sensorLoc,1)
        fh = mprfPlotHeadLayout(sensorLoc(ii,~isnan(sensorLoc(ii,:))),0,[]);
        
        if opt.saveFig
            aveName = regexp(dirPth.model.saveDataPth,'average', 'match');
            if isempty(aveName), aveName = dirPth.subjID; else aveName = aveName{1}; end
            if ~exist(fullfile(saveDir, type, 'sensorSelection'), 'dir'); mkdir(fullfile(saveDir, type, 'sensorSelection')); end
            print(gcf, fullfile(saveDir, type, 'sensorSelection', sprintf('fig2c_%s_vary%s_sensors%s_%s_%d', ...
                aveName, opt.vary.perturbOrigPRFs, opt.fNamePostFix, type, ii)), '-dpdf');
        end
        close(fh)
    end
else 
    fh = mprfPlotHeadLayout(sensorLoc,0,[]);
    if opt.saveFig
        if cellfind(strsplit(saveDir, '/'), 'average')
            fname = sprintf('fig2c_%s_vary%sSensors%s_%s', 'groupAverage', opt.fNamePostFix, opt.vary.perturbOrigPRFs, type);
        else
            fname = sprintf('fig2c_%s_vary%sSensors%s_%s', dirPth.subjID, opt.fNamePostFix, opt.vary.perturbOrigPRFs, type);
        end
        if ~exist(fullfile(saveDir, type,'sensorSelection'), 'dir'); mkdir(fullfile(saveDir,type, 'sensorSelection')); end
        print(gcf, fullfile(saveDir, type,'sensorSelection', fname), '-dpdf');
    end
    
    close(fh)
end

return
