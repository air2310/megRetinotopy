function plotSSVEFmesh(data, conditions, subjID, dirPth, opt)


%% Get dimensions of data, frequencies, and epochs of interest
[nSensors, nTimePoints, nEpochs] = size(data);
freqs = [opt.meg.flickerFreq-1,  opt.meg.flickerFreq, opt.meg.flickerFreq+1];

freqIdx = mprfFreq2Index(nTimePoints, freqs, opt.meg.fs);

epochsStimToPlot  = conditions<10;
epochsBlankToPlot = conditions==10;

%% INCOHERENT SPECTRUM
allAmps = abs(fft(data, [], 2))/size(data,2)*2;

% Get amplitudes of 10 Hz and divide by mean of surrounding (9 and 11 Hz)
% freq amplitudes, to get an estimate of SNR
amp10HzStim_incoh     = squeeze(nanmean(allAmps(:, freqIdx(2), epochsStimToPlot),3));
amp9to11HzStim_incoh  = nansum(nanmean(allAmps(:, [round(freqIdx(1)),  freqIdx(2), round(freqIdx(3))], epochsStimToPlot),3),2);
stimSNRToPlot.incoh   = amp10HzStim_incoh ./ amp9to11HzStim_incoh;

% Define regular stim and blank amplitude at 10 Hz for incoherent spectrum
stimDataToPlot.incoh    = amp10HzStim_incoh;
blankDataToPlot.incoh   = squeeze(nanmean(allAmps(:, freqIdx(2), epochsBlankToPlot),3));

%% COHERENT SPECTRUM
meanStimTs   = nanmean(data(:,:, epochsStimToPlot),3);
meanBlankTs  = nanmean(data(:,:, epochsBlankToPlot),3);

meanAmpsStim = abs(fft(meanStimTs, [], 2))/size(data,2)*2;
meanAmpsBlank = abs(fft(meanBlankTs, [], 2))/size(data,2)*2;

% freq amplitudes, to get an estimate of SNR
amp10HzStim_coh       = squeeze(meanAmpsStim(:, freqIdx(2)));
amp9to11HzStim_coh    = nansum(meanAmpsStim(:, [round(freqIdx(1)), freqIdx(2), round(freqIdx(3))]),2);
stimSNRToPlot.coh     =  amp10HzStim_coh ./ amp9to11HzStim_coh;

% Define regular stim and blank amplitude at 10 Hz for coherent spectrum
stimDataToPlot.coh = amp10HzStim_coh;
blankDataToPlot.coh =  squeeze(meanAmpsBlank(:,freqIdx(2)));


%% Plot meshes
fh1 = figure;
megPlotMap(stimSNRToPlot.coh,[min(stimSNRToPlot.coh) max(stimSNRToPlot.coh)],[],[],[],[],[],'interpmethod', 'nearest');
title('Steady state visually evoked field (10 Hz / sum 9-11 Hz); Coherent spectrum (max = 1)');

fh2 = figure;
megPlotMap(stimSNRToPlot.coh,[0 1],[],[],[],[],[],'interpmethod', 'nearest');
title('Steady state visually evoked field (10 Hz / sum 9-11 Hz); Coherent spectrum (max = 1)');

fh3 = figure;
megPlotMap(stimDataToPlot.coh-blankDataToPlot.coh,[0 max(stimDataToPlot.coh-blankDataToPlot.coh)],[],[],[],[],[],'interpmethod', 'nearest');
title('Steady state visually evoked field - blanks (10 Hz) Coherent spectrum');

fh4 = figure;
megPlotMap(stimSNRToPlot.incoh,[min(stimSNRToPlot.incoh) max(stimSNRToPlot.incoh)],[],[],[],[],[],'interpmethod', 'nearest');
title('Steady state visually evoked field (10 Hz / sum 9-11 Hz) Incoherent spectrum (max = 1)');

fh5 = figure;
megPlotMap(stimSNRToPlot.incoh,[0 1],[],[],[],[],[],'interpmethod', 'nearest');
title('Steady state visually evoked field (10 Hz / sum 9-11 Hz) Incoherent spectrum (max = 1)');


fh6 = figure;
megPlotMap(stimDataToPlot.incoh-blankDataToPlot.incoh,[0 max(stimDataToPlot.incoh-blankDataToPlot.incoh)],[],[],[],[],[],'interpmethod', 'nearest');
title('Steady state visually evoked field - blanks (10 Hz) Incoherent spectrum');


if opt.saveFig
    print(fh1, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_coh_SNR_maxCLim', subjID)))
    print(fh2, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_coh_SNR_0-1', subjID)))

    print(fh3, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_coh_diff', subjID)))
    print(fh4, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_incoh_SNR_maxCLim', subjID)))
    print(fh5, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_incoh_SNR_0-1', subjID)))
    print(fh6, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_incoh_diff', subjID)))
end

return