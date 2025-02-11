clear
clc
close all

%% Saving dirs
direct.dataroot = 'C:\Users\angel\Documents\MEGRetinotopicMappingData\';

subids = dir([direct.dataroot 'MEG\wlsub*']);
for sub = 1:length(subids)
    subjID = subids(sub).name;
    
    %% Setup paths 
    
    % Set options if not defined (see getOpts function for more options)
    if ~exist('opt', 'var') || isempty(opt)
        opt = getOpts('saveFig', true,'verbose', true, 'fullSizeMesh', true, ...
            'perturbOrigPRFs', false, 'addOffsetParam', false, ...
            'refitGainParam', false);
    end
    
    % Load paths with data files for this subject
    dirPth = loadPaths(subjID);
    
    % Go back to root
    cd(mprf_rootPath)
    
    fprintf('(%s): Starting analysis of subject %s, using %s\n', mfilename, subjID, regexprep(opt.subfolder,'/',' '));
    
    %% Load Preprocessed MEG data
    
    load(fullfile(dirPth.meg.processedDataPth, 'epoched_data_hp_preproc_denoised.mat'), 'data');
    load(fullfile(dirPth.meg.processedDataPth, 'megStimConditions.mat'), 'triggers');
    load(fullfile(dirPth.meg.processedDataPth, 'meg_stimulus.mat'), 'meg_stim');
    gainMtx    = loadGainMtx(subjID, dirPth, opt);
    
    % remove similar-named fields
    data       = data.data;
    stim       = meg_stim;
    conditions = triggers;
    
    % Get MEG data struct and fill with loaded data
    meg = struct();
    meg.data            = data;
    meg.stim            = stim;
    meg.stim.conditions = conditions;
    meg.gain            = gainMtx;
    
    % Save some memory
    clear gainMtx data conditions
    
    
    %% 2. Get MRI data path
    % This section will:
    % 2.1  Smoothing pRF params in voxel space (and recompute beta's)
    % 2.2  Export pRF params to FreeSurfer mesh + get Wang et al visual rois
    % 2.3  Export pRF to Brainstorm mesh (if opt.fullSizeMesh = false)
    
    % Get MRI data struct
    mri     = struct();
    
    % What pRF retinotopy maps are used on what size mesh?
    if opt.mri.useBensonMaps % Use benson maps or not
        [mri.data, mri.prfSurfPath] = loadBensonRetinotopyMaps(subjID, dirPth, opt);
    elseif opt.mri.useNYU3TAveMaps
        mri.prfSurfPath = loadNYU3TRetinotopyMaps(subjID, dirPth, opt);  
    elseif opt.fullSizeMesh % if using full gain matrix, we need pRF params on FreeSurfer surface
        mri.prfSurfPath = dirPth.fmri.saveDataPth_prfFS;
    else % if using downsampled gain matrix, we need pRF params on Brainstorm surface
        mri.prfSurfPath = dirPth.fmri.saveDataPth_prfBS;
    end
    
    
    %% 3. Forward model
    
    % 3.1 Predict response to MEG stimulus on mesh vertices (could be BS or FS)
    %       INPUTS  (1) path to pRF parameters on surface (string)
    %               (2) MEG stimulus (struct with x, y, images, etc)
    %       OUTPUTS (1) predicted responses on surface (epochs x vertices)
    
    predSurfResponse = mprf_MEGPredictionFromSurfaceWrapper(mri.prfSurfPath, meg.stim, dirPth, opt);
    results.predSurfResponse = predSurfResponse;
    
    
    %% 3.2 Predicted response for MEG stimulus at MEG sensor level (weighting
    %     predicted surface responses with gain matrix)
    %       INPUTS  (1) predicted surface responses (on BS or FS surface)
    %                                           (epochs x vertices)
    %               (2) gain matrix             (sensors x vertices)
    %
    %       OUTPUTS (1) predicted MEG responses (epochs x sensors)
    
    predMEGResponse = mprf_MEGPredictionSensorsWrapper(predSurfResponse, meg.gain, dirPth, opt);
    results.predMEGResponse = predMEGResponse;
        
    %% 3.3 Computing phase referenced amplitude from preprocessed MEG data
    % and predicted MEG responses from cortical surface
    %   	INPUTS  (1) preprocessed MEG data   (time x epochs x run x sensors)
    %               (2) predicted MEG responses (epochs x sensors)
    %       OUTPUTS (1) phase-referenced MEG time series
    %                                           (sensors x run group x epochs)
    % These data are rescaled according to how far the ssvep phase on each
    % trial was from a reference phase for which rescaling optomised regression
    % error between the MEG data and the predicted MEG response from fMRI. 
    %               (2) bestBetas               (1 x run group x sensors)
    %               (3) bestRefPhase            (1 x run group x sensors)
    %               (4) bestOffsets             (1 x run group x sensors)
    
    % [phaseRefMEGResponse, bestBetas, bestRefPhase, bestOffsets] = mprf_MEGPhaseReferenceDataWrapper(meg.data, predMEGResponse, 1, opt, dirPth);
    [phaseRefMEGResponse, bestBetas, bestRefPhase, bestOffsets] = mprf_MEGPhaseReferenceDataWrapper(meg.data, predMEGResponse, dirPth, opt);
    results.phaseRefMEGResponse = phaseRefMEGResponse; 
    results.bestBetas           = bestBetas;
    results.bestRefPhase        = bestRefPhase;
    results.bestOffsets         = bestOffsets;

    %% Stim properties
    %numepochs = 140
    % nruns = 19
    %numverticels = 291974
    %numsensors = 157
    %numtimepoints = 1100
    % opt.meg.fs=        1000
    %megflickerfreq = 10hz
    %freqidx=12
    % size(meg.dta) =
    % timepoints   numepochs    nruns      numsensors
    % 1100         140          19         157
    

    % stim = meg.stim;
    numepochs = size(meg.data,2);
    stimuse = meg.stim.resizedIm;
    epochsuse = not(meg.stim.conditions.stimConditions(1:140)==20);
    % for epoch = 1:numepochs
    %     idx = find(stim.im(:,epoch));
    %     plot(stim.X(idx), stim.Y(idx), 'x')
    % end
    
    %% Plot resposes by epoch. 
    
    figure;
    datplot = mean(phaseRefMEGResponse(epochsuse,:,:),2);
    datplot(isnan(datplot)) = 0;
    index = 6:2:25;%61:2:90;
    clims = [min(min(datplot(index,:))), max(max(datplot(index,:)))];
    count = 0;
    for ee = index
        count = count + 1;
        subplot(2,5,count)    
        megPlotMap(squeeze(datplot(ee,:)),clims,[],[],[],[],[],'interpmethod', 'nearest');
    end
    title(subjID)
    
    figure;
    count = 0;
    for ee = index
        count = count + 1;
        subplot(2,5,count)    
        imagesc(stimuse(:,:,ee))
    end
    
    
    % what about the raw data?
    % figure;
    % megPlotMap(squeeze( nanmean(nanmean(meg.data(:,10,:,:),3),1)),[],[],[],[],[],[],'interpmethod', 'nearest');
    
    %% Save for python
     conds = meg.stim.conditions.stimConditions(1:140);

     save([direct.dataroot 'PythonData\' subjID  '.mat'], 'phaseRefMEGResponse')
     save([direct.dataroot 'PythonData\' subjID  'stim.mat'], 'stimuse')
     save([direct.dataroot 'PythonData\' subjID  'cond.mat'], 'conds')

end

%% Create Grand average images
allresponses = NaN(140,157,30);
for sub = 1:length(subids)
    subjID = subids(sub).name;

    dat = load([direct.dataroot 'PythonData\' subjID  '.mat']);
    allresponses(:,:,sub) = squeeze(nanmean(dat.phaseRefMEGResponse,2));

end

% plot

h=figure;
datplot = nanmean(allresponses,3);
datplot(isnan(datplot)) = 0;
index = 6:2:25;%61:2:90;
clims = [min(min(datplot(index,:))), max(max(datplot(index,:)))];
count = 0;
for ee = index
    count = count + 1;
    subplot(2,5,count)    
    megPlotMap(squeeze(datplot(ee,:)),clims,[],[],[],[],[],'interpmethod', 'nearest');
    title(['Position ', num2str(ee)])
end
legend('off')
suptitle('Subject Mean MEG Response')
saveas(h, [direct.dataroot 'PythonData\MEGSensorResponse.png'])