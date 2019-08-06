function plotSSVEFmesh(data, conditions, subjID, dirPth, opt)


%% Get dimensions of data, frequencies, and epochs of interest
[nSensors, nTimePoints, nEpochs] = size(data);
freqs = [opt.flickerFreq-1,  opt.flickerFreq, opt.flickerFreq+1];

freqIdx = mprfFreq2Index(nTimePoints, freqs, opt.fs);

epochsStimToPlot  = conditions<10;
epochsBlankToPlot = conditions==10;

%% INCOHERENT SPECTRUM
allAmps = abs(fft(data, [], 2))/size(data,2)*2;

% Get amplitudes of 10 Hz and divide by mean of surrounding (9 and 11 Hz)
% freq amplitudes, to get an estimate of SNR
amp10HzStim_incoh     = squeeze(nanmean(allAmps(:, freqIdx(2), epochsStimToPlot),3));
amp9and11HzStim_incoh = squeeze(nanmean(allAmps(:, [round(freqIdx(1)), round(freqIdx(3))], epochsStimToPlot),3));
stimSNRToPlot.incoh   = amp10HzStim_incoh ./ nanmean(amp9and11HzStim_incoh,2);

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
amp9and11HzStim_coh   = squeeze(nanmean(meanAmpsStim(:, [round(freqIdx(1)), round(freqIdx(3))]),2));
stimSNRToPlot.coh     =  amp10HzStim_coh ./ nanmean(amp9and11HzStim_coh,1);

% Define regular stim and blank amplitude at 10 Hz for coherent spectrum
stimDataToPlot.coh = amp10HzStim_coh;
blankDataToPlot.coh =  squeeze(meanAmpsBlank(:,freqIdx(2)));


%% Plot meshes
fh1 = figure;
megPlotMap(stimSNRToPlot.coh,[0 max(stimSNRToPlot.coh)],[],[],[],[],[],'interpmethod', 'nearest');
title('Steady state visually evoked field (10 Hz ./ +-1 Hz); Coherent spectrum');

fh2 = figure;
megPlotMap(stimDataToPlot.coh-blankDataToPlot.coh,[0 max(stimDataToPlot.coh-blankDataToPlot.coh)],[],[],[],[],[],'interpmethod', 'nearest');
title('Steady state visually evoked field - blanks (10 Hz) Coherent spectrum');

fh3 = figure;
megPlotMap(stimSNRToPlot.incoh,[1 max(stimSNRToPlot.incoh)],[],[],[],[],[],'interpmethod', 'nearest');
title('Steady state visually evoked field (10 Hz ./ +-1 Hz) Incoherent spectrum');

fh4 = figure;
megPlotMap(stimDataToPlot.incoh-blankDataToPlot.incoh,[0 max(stimDataToPlot.incoh-blankDataToPlot.incoh)],[],[],[],[],[],'interpmethod', 'nearest');
title('Steady state visually evoked field - blanks (10 Hz) Incoherent spectrum');


if opt.saveFig
    print(fh1, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_coh_SNR', subjID)))
    print(fh2, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_coh_diff', subjID)))
    print(fh3, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_incoh_SNR', subjID)))
    print(fh4, '-dpng', fullfile(dirPth.meg.saveFigPth,sprintf('%s_SSVEFMESH_postDenoise_incoh_diff', subjID)))
end

return