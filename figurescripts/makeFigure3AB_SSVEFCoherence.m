function [] = makeFigure3AB_SSVEFCoherence()
% Function to plot individual subjects for Figure 3A and B from manuscript,
% plotting 10 Hz SSVEF coherence (10Hz/(9-11Hz).

subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};


fH1 = figure(1); clf; set(gcf,'Position',[0,300,500,500]); set(fH1, 'Name', 'SSVEF coherent spectrum' , 'NumberTitle', 'off');
fH2 = figure(2); clf; set(gcf,'Position',[0,300,500,500]); set(fH2, 'Name', 'SSVEF incoherent spectrum' , 'NumberTitle', 'off');
fH3 = figure(3); clf; set(gcf, 'Position',[0,300,500,500]); set(fH3, 'Name', 'Spectrum of one MEG sensor' , 'NumberTitle', 'off');

saveSubDir = 'Figure3_SSVEF_coherence';
saveDir = fullfile(mprf_rootPath,'data','Retinotopy','Quality_check','average','finalfig',saveSubDir);

interpMethod = 'v4'; % or if not interpolated, use 'nearest'

if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

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
    meanBlankTs  = nanmean(dataReshaped(:,:, epochsBlankToPlot),3);
    
    meanAmpsStim = abs(fft(meanStimTs, [], 2))/size(dataReshaped,2)*2;
    meanAmpsBlank = abs(fft(meanBlankTs, [], 2))/size(dataReshaped,2)*2;
    
    % freq amplitudes, to get an estimate of SNR
    amp10HzStim_coh       = squeeze(meanAmpsStim(:, freqIdx(2)));
    amp9to11HzStim_coh    = nansum(meanAmpsStim(:, [round(freqIdx(1)), freqIdx(2), round(freqIdx(3))]),2);
    stimSNRToPlot_coh(s,:)  =  amp10HzStim_coh ./ amp9to11HzStim_coh;
    
    %% Plot spectra
    figure(fH3); clf;
    
    % Define MEG sensor to plot
    if s == 1
        sensorIdx = 13;
    else
        sensorIdx = 1;
    end
    
    % Define color to plot full conditions
    colors    = [0 0 0; 126 126 126]./255;
    
    % Spectrum of example channel
    t = (0:size(meanAmpsStim,2)-1)./opt.meg.fs;
    f = (0:length(t)-1)/max(t);    
    xl = [5 60];    
    fok = f;
    fok(f==60) = [];
    xt = [10:10:60];
    ymax = 0.05+(max(meanAmpsStim(sensorIdx,10:50)).*10^14);
    yt = [0:.1:ymax]; % .*(10^-14) => fempto tesla
    yl = [yt(1), yt(end)]; % .*(10^-14) => fempto tesla
    
    % plot mean
    plot(f, meanAmpsStim(sensorIdx,:),  '-',  'Color', colors(1,:), 'LineWidth', 2); hold on;
    plot(f, meanAmpsBlank(sensorIdx,:),  '-',  'Color', colors(2,:), 'LineWidth', 2);
    
    % format x and y axes
    set(gca, 'XLim', xl, 'XTick', xt, 'XScale', 'linear', 'YScale','linear');
    set(gca,'ytick',yt.*(10^-14),'yticklabel', sprintfc('%1.1f',yt), 'ylim',yl.*(10^-14), 'TickDir', 'out', 'FontSize', 12);
    box off;
    
    % label figure, add stimulus harmonic lines, and make it look nice
    xlabel('Frequency (Hz)'); ylabel('Amplitude (fT/Hz)');
    title(sprintf('Coherent spectrum sensor %d', sensorIdx));
    yl2 = get(gca, 'YLim');
    for ii = xt, plot([ii ii], yl2, '-', 'Color', [100 100 100]./255); end
    legend({'Stimulus (10 Hz contrast reversals)', 'blank'}, 'Location', 'NorthEast'); legend boxoff
    set(gca, 'FontSize', 12)
    if opt.saveFig
        figurewrite(fullfile(saveDir, sprintf('Figure3_SpectrumOneChannel_S%d_usingCoherentSpectrum_%s', s, opt.fNamePostFix)), [],0,'.',1);
        print(fullfile(saveDir, sprintf('Figure3_SpectrumOneChannel_S%d_usingCoherentSpectrum_%s', s, opt.fNamePostFix)), '-dpng');
    end
    
    %% Plot meshes
    figure(fH1);
    subplot(2,5,s)
    megPlotMap(stimSNRToPlot_coh(s,:),[0 .9], fH1, 'parula',[],[],[],'interpmethod', interpMethod);
    title(sprintf('S%d',s));
    c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
    pos = c.Position; set(c, 'Position', [pos(1)+0.09 pos(2)+0.06, pos(3)/1.5, pos(4)/1.5])
    %         figurewrite(fullfile(saveDir, sprintf('Figure3_S%d_SSVEFcoherence_usingCoherentSpectrum_%s_%s', s, opt.fNamePostFix,interpMethod)),[],0,'.',1);

    figure(fH2); 
    subplot(2,5,s) 
    megPlotMap(stimSNRToPlot_incoh(s,:),[.3 .4], fH1, 'parula',[],[],[],'interpmethod', interpMethod);
    title(sprintf('S%d',s));
    c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
    pos = c.Position; set(c, 'Position', [pos(1)+0.09 pos(2)+0.06, pos(3)/1.5, pos(4)/1.5])
    %         figurewrite(fullfile(saveDir, sprintf('Figure3_S%d_SSVEFcoherence_usingIncoherentSpectrum_%s_%s', s, opt.fNamePostFix,interpMethod)),[],0,'.',1);
end
   
% Save all subjects in one figure
if opt.saveFig
    figure(fH1)
    print(fH1, fullfile(saveDir, sprintf('Figure3_SSVEFcoherence_usingCoherentSpectrum_%s_%s',opt.fNamePostFix,interpMethod)), '-dpng');
    print(fH1, fullfile(saveDir, sprintf('Figure3_SSVEFcoherence_usingCoherentSpectrum_%s_%s',opt.fNamePostFix,interpMethod)), '-dpdf');
    fprintf('\n saving figure 3: All subject SSVEF_coherence in %s',saveDir);

    figure(fH2)
    print(fH2, fullfile(saveDir, sprintf('Figure3_SSVEFcoherence_usingIncoherentSpectrum_%s_%s', opt.fNamePostFix,interpMethod)), '-dpng');
    print(fH2, fullfile(saveDir, sprintf('Figure3_SSVEFcoherence_usingIncoherentSpectrum_%s_%s', opt.fNamePostFix,interpMethod)), '-dpdf');
    fprintf('\n saving figure 3: All subjects SSVEF_coherence in %s',saveDir);
end



fH4 = figure();
mn_stimSNRToPlot_coh = nanmean(stimSNRToPlot_coh,1);
megPlotMap(mn_stimSNRToPlot_coh,[0 .9], fH2, 'parula',[],[],[],'interpmethod', interpMethod);
title(sprintf('Group Average (N=%d) SSVEF coherence: Incoherent spectrum',length(subjects)));
c = colorbar; c.TickDirection = 'out'; c.Box = 'off';

fH5 = figure();
mn_stimSNRToPlot_incoh = nanmean(stimSNRToPlot_incoh,1);
megPlotMap(mn_stimSNRToPlot_incoh,[.3 .4], fH2, 'parula',[],[],[],'interpmethod', interpMethod);
title(sprintf('Group Average (N=%d) SSVEF coherence: Incoherent spectrum',length(subjects)));
c = colorbar; c.TickDirection = 'out'; c.Box = 'off';

if opt.saveFig
    figure(fH4)
    print(fH4, fullfile(saveDir, sprintf('Figure3_GroupAverage_SSVEFcoherence_usingCoherentSpectrum_%s_%s', opt.fNamePostFix)), '-dpng');
    figurewrite(fullfile(saveDir, sprintf('Figure3_GroupAverage_SSVEFcoherence_usingCoherentSpectrum_%s_%s', opt.fNamePostFix)),[],0,'.',1);
    fprintf('\n saving figure XX: Group Average SSVEF_coherence in %s',saveDir);
    
    figure(fH5)
    print(fH5, fullfile(saveDir, sprintf('Figure3_GroupAverage_SSVEFcoherence_usingIncoherentSpectrum_%s_%s', opt.fNamePostFix)), '-dpng');
    figurewrite(fullfile(saveDir, sprintf('Figure3_GroupAverage_SSVEFcoherence_usingIncoherentSpectrum_%s_%s', opt.fNamePostFix)),[],0,'.',1);
    fprintf('\n saving figure XX: Group Average SSVEF_coherence in %s',saveDir);
end