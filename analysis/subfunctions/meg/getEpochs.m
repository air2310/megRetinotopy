function epochedData = getEpochs(data, triggers, epochLength, flickerFreq, fs)


% Get epoch samples
numTimepoints  = round((epochLength(2)-epochLength(1))*fs);
skipTimepoints = round(epochLength(1)*fs); % Question: skip 150 ms of the first bar location after a blank or for every bar position?
numChannels    = size(data,1);
epochSamples   = [skipTimepoints, numTimepoints+skipTimepoints]; %epoch length in samples


% Not USED: Skip every 10 triggers, to the triggers that define the start of every run
% endOfRun = floor(length(triggers.timing)/flickerFreq : length(triggers.timing)/flickerFreq : length(triggers.timing));
% startOfRun = [1, 1+endOfRun(1:end-1)];


% Get indices for epoching
inds = bsxfun(@plus,triggers.timing,(epochSamples(1):epochSamples(2)-1));

% Epoch data and reshape into time x epochs x channels
epochedData   = data(:,inds); 
epochedData   = reshape(epochedData, numChannels, size(inds,1),size(inds,2));
epochedData   = permute(epochedData, [3 2 1]); % time points x epochs x channels

end

















