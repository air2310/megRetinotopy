function ab = mprfDenoiseEvalFun(data, F,s_rate)


tmp =  (size(data,2)./(s_rate./F)) + 1;
keep_idx = round(tmp(1)) : round(tmp(2));

% Fourier transform the data, use only the relevant frequencies
spec = fft(data,[],2);
spec_amp = abs(spec(:,keep_idx,:))/ size(data,2)*2; 

% take the mean across frequencies in log space 
ab = squeeze(exp(nanmean(log(spec_amp.^2),2)))';






end








