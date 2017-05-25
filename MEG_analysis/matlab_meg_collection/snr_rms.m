function [snrvector, RMS] = snr_rms(data, pre, smpr)
pre = pre * smpr;
RMS = sqrt(mean((data(:, :).^2), 2));
snrvector = RMS ./ mean(RMS(1:pre));