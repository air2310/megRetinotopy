function sl = mprfDenoiseEvalFun(data, F,s_rate)

tmp =  (size(data,2)./(s_rate./F)) + 1;
keepIdx = round(tmp(1)) : round(tmp(2));

% Fourier transform the data, use only the relevant frequencies
spectralData = fft(data,[],2);
spectralAmplitude = abs(spectralData(:,keepIdx,:))/ size(data,2)*2; 

% take the mean across frequencies in log space 
sl = squeeze(exp(nanmean(log(spectralAmplitude.^2),2)))';

end








