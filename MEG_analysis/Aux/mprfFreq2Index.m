function index = mprfFreq2Index(n_samples, F, s_rate)

index = (n_samples./(s_rate./F)) + 1;


end