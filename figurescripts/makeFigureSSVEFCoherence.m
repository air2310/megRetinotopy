function [] = makeFigureSSVEFCoherence()

subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
            'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

          
fH1 = figure; set(gcf,'Position',[1000, 651, 1500, 687]);
fH2 = figure; set(gcf,'Position',[1000, 651, 1500, 687]);


for s = 1:length(subjects)

    % Get subject ID, options and paths
    subjID = subjects{s};
    opt    = getOpts('verbose', true, 'saveFig',true, 'headmodel', 'OS','fullSizeMesh',true);
    dirPth = loadPaths(subjID);

    % Load data files for this subject
    load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'epoched_data_hp_preproc_denoised.mat'), 'data');
    load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'megStimConditions.mat'), 'triggers');
    
    % Remove similar named field
    data = data.data;
    
    % Permute data
    dataPermuted = permute(data, [4 1 2 3]);
    
    % Get array size of MEG data 
    [nSensors, nTimePoints, nEpochs, nRuns] = size(dataPermuted);

    % Reshape data
    dataReshaped =reshape(dataPermuted, [nSensors, nTimePoints, nEpochs*nRuns]);

    % Get SSVEF frequency of stim and epochs with stim
    freqs = [opt.meg.flickerFreq-1,  opt.meg.flickerFreq, opt.meg.flickerFreq+1];
    freqIdx = mprfFreq2Index(nTimePoints, freqs, opt.meg.fs);

    epochsStimToPlot  = triggers.stimConditions<10;
    epochsBlankToPlot = triggers.stimConditions==10;

    %% Compute INCOHERENT SPECTRUM
    allAmps = abs(fft(dataReshaped, [], 2))/size(dataReshaped,2)*2;

    % Get amplitudes of 10 Hz and divide by mean of surrounding (9 and 11 Hz)
    % freq amplitudes, to get an estimate of SNR
    amp10HzStim_incoh     = squeeze(nanmean(allAmps(:, freqIdx(2), epochsStimToPlot),3));
    amp9to11HzStim_incoh  = nansum(nanmean(allAmps(:, [round(freqIdx(1)),  freqIdx(2), round(freqIdx(3))], epochsStimToPlot),3),2);
    stimSNRToPlot_incoh(s,:)   = amp10HzStim_incoh ./ amp9to11HzStim_incoh;

    %% Compute COHERENT SPECTRUM
    meanStimTs   = nanmean(dataReshaped(:,:, epochsStimToPlot),3);
    meanAmpsStim = abs(fft(meanStimTs, [], 2))/size(dataReshaped,2)*2;

    % freq amplitudes, to get an estimate of SNR
    amp10HzStim_coh       = squeeze(meanAmpsStim(:, freqIdx(2)));
    amp9to11HzStim_coh    = nansum(meanAmpsStim(:, [round(freqIdx(1)), freqIdx(2), round(freqIdx(3))]),2);
    stimSNRToPlot_coh(s,:)  =  amp10HzStim_coh ./ amp9to11HzStim_coh;



    %% Plot meshes
    figure(fH1)
    subplot(2,5,s)
    megPlotMap(stimSNRToPlot_coh(s,:),[0 1], fH1, 'parula',[],[],[],'interpmethod', 'nearest');
    title(sprintf('S%d',s));
    c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
    pos = c.Position; set(c, 'Position', [pos(1)+0.02 pos(2)+0.07, pos(3)/2, pos(4)/2])
    
    figure(fH2)
    subplot(2,5,s)
    megPlotMap(stimSNRToPlot_incoh(s,:),[0.3 .6], fH1, 'parula',[],[],[],'interpmethod', 'nearest');
    title(sprintf('S%d',s));
    c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
    pos = c.Position; set(c, 'Position', [pos(1)+0.02 pos(2)+0.07, pos(3)/2, pos(4)/2])
    
end

fH3 = figure();
mn_stimSNRToPlot_coh = nanmean(stimSNRToPlot_coh,1);
megPlotMap(mn_stimSNRToPlot_coh,[0 1], fH2, 'parula',[],[],[],'interpmethod', 'nearest');
title(sprintf('Group Average (N=%d) SSVEF coherence: Incoherent spectrum',length(subjects)));
c = colorbar; c.TickDirection = 'out'; c.Box = 'off';

fH4 = figure();
mn_stimSNRToPlot_incoh = nanmean(stimSNRToPlot_coh,1);
megPlotMap(mn_stimSNRToPlot_coh,[0 1], fH2, 'parula',[],[],[],'interpmethod', 'nearest');
title(sprintf('Group Average (N=%d) SSVEF coherence: Incoherent spectrum',length(subjects)));
c = colorbar; c.TickDirection = 'out'; c.Box = 'off';

if opt.saveFig
    saveSubDir = 'figureXX_SSVEF_coherence';
    saveDir = fullfile(mprf_rootPath,'data','Retinotopy','Quality_check','average','finalfig',saveSubDir);
   
    print(fH1, fullfile(saveDir, sprintf('figXX_IndividualSubjects_SSVEFcoherence_usingCoherentSpectrum_%s', opt.fNamePostFix)), '-dpng');
    print(fH1, fullfile(saveDir, sprintf('figXX_IndividualSubjects_SSVEFcoherence_usingCoherentSpectrum_%s', opt.fNamePostFix)), '-depsc');
    fprintf('\n saving figure XX: Individual Subject SSVEF_coherence in %s',saveDir);
   
    print(fH2, fullfile(saveDir, sprintf('figXX_IndividualSubjects_SSVEFcoherence_usingIncoherentSpectrum_%s', opt.fNamePostFix)), '-dpng');
    print(fH2, fullfile(saveDir, sprintf('figXX_IndividualSubjects_SSVEFcoherence_usingIncoherentSpectrum_%s', opt.fNamePostFix)), '-depsc');
    fprintf('\n saving figure XX: Individual Subject SSVEF_coherence in %s',saveDir);
   
    
    print(fH3, fullfile(saveDir, sprintf('figXX_GroupAverage_SSVEFcoherence_usingCoherentSpectrum_%s', opt.fNamePostFix)), '-dpng');
    print(fH3, fullfile(saveDir, sprintf('figXX_GroupAverage_SSVEFcoherence_usingCoherentSpectrum_%s', opt.fNamePostFix)), '-depsc');
    fprintf('\n saving figure XX: Group Average SSVEF_coherence in %s',saveDir);
    
    print(fH4, fullfile(saveDir, sprintf('figXX_GroupAverage_SSVEFcoherence_usingIncoherentSpectrum_%s', opt.fNamePostFix)), '-dpng');
    print(fH4, fullfile(saveDir, sprintf('figXX_GroupAverage_SSVEFcoherence_usingIncoherentSpectrum_%s', opt.fNamePostFix)), '-depsc');
    fprintf('\n saving figure XX: Group Average SSVEF_coherence in %s',saveDir);
end