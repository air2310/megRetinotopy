function [] = makeFigure3AB_SSVEFCoherence(subjectToPlot,plotAverage,plotSupplementalFig)
% Function to plot individual subjects for Figure 3A & B from manuscript,
% plotting 10 Hz SSVEF coherence (10Hz/(9-11Hz).
%
% INPUTS:
% subjectsToPlot      :  Vector with integers to select subjects to plot
%                         should be between 1-10.
% plotAverage         :  Boolean to plot additonal figure with group average.
% plotSupplementalFig :  Boolean to plot supplementary figure with all
%                         subjects.
%
%

subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

figSize = [0,300,500,500];
if plotAverage
    subjectToPlot = 1:length(subjects);
    saveSubDir = 'Figure3AB_SSVEF_coherence';
elseif plotSupplementalFig
    subjectToPlot = 1:length(subjects);
    saveSubDir = 'SupplFigure1_SSVEF_coherence';
    figSize = get(0, 'screensize');
else
   saveSubDir = 'Figure3AB_SSVEF_coherence';
   fH3 = figure(3); clf; set(gcf, 'Position',figSize); set(fH3, 'Name', 'Spectrum of one MEG sensor' , 'NumberTitle', 'off');
end

fH1 = figure(1); clf; set(gcf,'Position',figSize); set(fH1, 'Name', 'SSVEF coherent spectrum' , 'NumberTitle', 'off');

% Plotting params
interpMethod = 'v4'; % or if not interpolated, use 'nearest'

if length(subjectToPlot) > 5
    nrows = 2;
    ncols = 5;
else
    nrows = 1;
    ncols = length(subjectToPlot);
end

for s = subjectToPlot
    
    % Get subject ID, options and paths
    subjID  = subjects{s};
    opt     = getOpts('verbose', true, 'saveFig',true, 'headmodel', 'OS','fullSizeMesh',true);
    dirPth  = loadPaths(subjID);
    if plotSupplementalFig
        saveDir = fullfile(dirPth.finalFig.savePthAverage,saveSubDir);
    else
        saveDir = fullfile(dirPth.finalFig.savePth,saveSubDir);
    end
    if ~exist(saveDir,'dir')
        mkdir(saveDir);
    end
    
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
    
%     %% Compute INCOHERENT SPECTRUM
%     allAmps = abs(fft(dataReshaped, [], 2))/size(dataReshaped,2)*2;
%     
%     % Get amplitudes of 10 Hz and divide by mean of surrounding (9 and 11 Hz)
%     % freq amplitudes, to get an estimate of SNR
%     amp10HzStim_incoh     = squeeze(nanmean(allAmps(:, freqIdx(2), epochsStimToPlot),3));
%     amp9to11HzStim_incoh  = nansum(nanmean(allAmps(:, [round(freqIdx(1)),  freqIdx(2), round(freqIdx(3))], epochsStimToPlot),3),2);
%     stimSNRToPlot_incoh(s,:)   = amp10HzStim_incoh ./ amp9to11HzStim_incoh;
%     
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
    if ~(plotSupplementalFig || plotAverage)
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
            figurewrite(fullfile(saveDir, sprintf('Figure3A_SpectrumOneChannel_S%d_usingCoherentSpectrum_%s', s, opt.fNamePostFix)), [],0,'.',1);
            print(fullfile(saveDir, sprintf('Figure3A_SpectrumOneChannel_S%d_usingCoherentSpectrum_%s', s, opt.fNamePostFix)), '-dpng');
        end
    end
    
    %% Plot meshes
    figure(fH1);
    subplot(nrows,ncols,s)
    megPlotMap(stimSNRToPlot_coh(s,:),[0 .9], fH1, 'parula',[],[],[],'interpmethod', interpMethod);
    title(sprintf('S%d',s));
    if (s==ncols) || ((2*s)==(2*ncols))
        c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
        pos = c.Position; set(c, 'Position', [pos(1)+0.09 pos(2)+0.06, pos(3)/1.5, pos(4)/1.5])
    else
        colorbar off
    end
end

% Save figure if requested
if opt.saveFig
    figure(fH1)
    print(fH1, fullfile(saveDir, sprintf('Figure3B_SSVEFcoherence_usingCoherentSpectrum_%s_%s',opt.fNamePostFix,interpMethod)), '-dpng');
    print(fH1, fullfile(saveDir, sprintf('Figure3B_SSVEFcoherence_usingCoherentSpectrum_%s_%s',opt.fNamePostFix,interpMethod)), '-dpdf');
    fprintf('(%s) Saving Figure 3B: Subject SSVEF_coherence in %s\n',mfilename,saveDir);
end


if plotAverage
    fH4 = figure(4); clf; set(gcf, 'Position',[0,300,500,500]); set(fH4, 'Name', 'SSVEF Coherence Group Average' , 'NumberTitle', 'off');
    
    mn_stimSNRToPlot_coh = nanmean(stimSNRToPlot_coh,1);
    megPlotMap(mn_stimSNRToPlot_coh,[0 .9], fH4, 'parula',[],[],[],'interpmethod', interpMethod);
    title(sprintf('Group Average (N=%d) SSVEF coherence: using coherent spectrum',length(subjects)));
    c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
    
    if opt.saveFig
        saveDir = fullfile(dirPth.finalFig.savePthAverage,saveSubDir);
        if ~exist(saveDir, 'dir')
            mkdir(saveDir)
        end
        figure(fH4)
        print(fH4, fullfile(saveDir, sprintf('Figure3B_GroupAverage_SSVEFcoherence_usingCoherentSpectrum_%s_%s', opt.fNamePostFix)), '-dpng');
        figurewrite(fullfile(saveDir, sprintf('Figure3B_GroupAverage_SSVEFcoherence_usingCoherentSpectrum_%s_%s', opt.fNamePostFix)),[],0,'.',1);
        fprintf('(%s) Saving Figure 3B: Group Average SSVEF coherence in %s \n',mfilename,saveDir);
    end
end