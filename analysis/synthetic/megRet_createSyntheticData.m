function megRet_createSyntheticData()
% This is a function to create a synthetic MEG dataset to test the analysis
% pipeline of the megRetinotopy project. 

%% 0. Define parameters
numTimePoints = 1100;       % time points per epoch
fs = 1000;                  % Sample rate (Hz)
t  = (1:numTimePoints)/fs;  % time vector in ms

numChan   = 157; % number of meg sensors
numRuns   = 19;  % number of runs
numEpochs = 140; % number of 1300-ms epochs

A_SL    = 10; % amplitude of stimulus-locked signal
A_BB    = 0; % amplitude of broadband response
A_noise = 1; % amplitude of noise response

% save data?
saveData = false;
plotFigs = true;

%% Load Yokogawa MEG header for channel positions
load('meg160_example_hdr');
xyz = hdr.grad.chanpos;

% Find the 75 channels in front (highest x value): these will be noise pool
[~, sortY] = sort(xyz(:,1), 'descend');
sensors.upperhalf = sortY(1:75);
sensors.lowerhalf = sortY(76:end); 

% Among the 82 non-noise pool, find the 41 channels that are most left and
% the 41 that are most right
[~, sortX] = sort(xyz(sensors.lowerhalf,2), 'descend');
sensors.lowerright  = sensors.lowerhalf(sortX(1:41));
sensors.lowerleft   = sensors.lowerhalf(sortX(42:end));

% Create map
map = zeros(1, 157);
map(sensors.upperhalf) = 1;
map(sensors.lowerleft) = 2;
map(sensors.lowerright) = 3;

% Plot it
megPlotMap(map)

%% 1. Generate Stim-Locked
stimLockedSignal = A_SL.*sin(t*2*pi*10); % 10 Hz sine wave

%% 2. Generate Pink noise
f       = -numTimePoints/2:numTimePoints/2-1;
alpha   = 1;
cpfov   = ifftshift(f);
amp     = 1./(abs(cpfov).^alpha);
amp(1)  = 0; % remove DC

%% 3. load conditions

tmp = load(fullfile(mprfRootPath, 'Synthetic', 'megStimConditions'), 'triggers');
conditions = tmp.triggers.stimConditions; %10=blink, 20=blank, 1,3,4,6,7 = 0, 90, 135, 225, 270 deg orientation
% we'll only use three conditions:
% 0 deg - lower left sensors
% 90 deg - lower right sensors,
% 225 deg - both lower left and right sensors
% 135, 270 will be noise controls)

data = NaN(numEpochs, numRuns, numChan, numTimePoints);


for ii = 1:numRuns
    fprintf('.');
    for jj = 1:numEpochs
       
        thisCondition = conditions(jj*ii);
        switch thisCondition
            case 1
                chanSignal = sensors.lowerleft;
                chanNoise  = setdiff(1:numChan, chanSignal);
                
            case 3
                 chanSignal = sensors.lowerright;
                 chanNoise  = setdiff(1:numChan, chanSignal);
               
            case 6
                chanSignal = sensors.lowerhalf;
                chanNoise  = setdiff(1:numChan, chanSignal);
                
            case {4,7, 10, 20}
                chanSignal = [];
                chanNoise  = 1:numChan;

        end
        
        % Simulate noise
        backgroundNoise = real(ifft(fft(randn(numTimePoints,length(chanNoise)),[],1) .* repmat(amp',[1 length(chanNoise)]),[],1));
        data(jj, ii, chanNoise, :)  = A_noise.*backgroundNoise';
        
        % Simulate signal + noise
        if ~isempty(chanSignal)
            sl    = repmat(stimLockedSignal',[1 length(chanSignal)]);
            bb    = real(ifft(fft(randn(numTimePoints,length(chanSignal)),[],1) .* repmat(amp',[1 length(chanSignal)]),[],1));
            noise = real(ifft(fft(randn(numTimePoints,length(chanSignal)),[],1) .* repmat(amp',[1 length(chanSignal)]),[],1));

            bb      = bb' * A_BB;
            noise   = noise' * A_noise;
            
            sl = sl';
            combinedData = sl+bb+noise;
            combinedData = bsxfun(@rdivide, combinedData, max(combinedData,[],2));
            data(jj, ii, chanSignal, :)  = combinedData;
        end
       
        
        
    end
end
data = permute(data, [4, 1, 3, 2]);

if saveData
    save(fullfile(mprfRootPath, 'Synthetic', 'meg_data'), 'data', '-v7.3')
end

if plotFigs
    conditions = reshape(conditions, [140,19]);
    % Plot mean epoch time series of a single condition
    figure; plot(squeeze(data(:,conditions(:,1)==1,1,1)));
    figure; plot(squeeze(data(:,5,1,1)));

    amps = abs(fft(data(:,conditions(:,1)==1,:,1)));
    figure; plot(0:1099, squeeze(mean(amps,2))); set(gca, 'XLim', [1 150], 'YScale', 'log')
    ylabel('Amplitude (AU)'); xlabel('Frequency (Hz)')
    figure; megPlotMap(squeeze(mean(amps(11,:,:),2)))
    
end
end
