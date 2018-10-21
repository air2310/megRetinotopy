function epochedData = getEpochs(data, triggers, epochLength, flickerFreq, fs)

numTimepoints  = epochLength(2)-epochLength(1);
skipTimepoints = epochLength(1); % Question: skip 150 ms of the first bar location after a blank or for every bar position?
numChannels    = size(data,1);


% Skip every 10 triggers, to the triggers that define the start of every
% run
endOfRun = floor(length(triggers.timing)/flickerFreq : length(triggers.timing)/flickerFreq : length(triggers.timing));
startOfRun = [1, 1+endOfRun(1:end-1)];

epochSamples = round(numTimepoints * fs); %epoch length in samples




inds = bsxfun(@plus,triggers.timing,(epochLength(1):epochLength(2)-1));



nEpochs = [];

% Preallocate space for epoched data array
epochedData = nan(numTimepoints, nEpochs, numChannels);

epoch_samples = round(epoch_time * fs); %epoch length in samples

inds = bsxfun(@plus,onset_times,(epoch_samples(1):epoch_samples(2)-1));

ts   = raw_ts(inds,:); 
ts   = reshape(ts,size(inds,1),size(inds,2),size(raw_ts,2));
ts   = permute(ts, [2 1 3]);



cur_step = 1;
for n = 1:length(startOfRun)
    
    cur_per = mod(n,n_periods);
    cur_rep = ceil(n/n_periods);
    
    if cur_per == 0
        cur_per = n_periods;
    end
    
    
    if mod(n,fb(cur_step)) == 0
        n_perc = n / length(triggers) * 100;
        fprintf('%0.2f%%\n',n_perc)
        cur_step = cur_step + 1;
    end
    
    if n <= n_periods
        epoch_id(n) = triggers(n,2);
    end
    
    start_idx = triggers(n,1)+skip_n_first_data_points;
    
    
    if n_data_points_per_epoch == -1
        if cur_per == n_periods % At the end of a period, take 1300 ms
            end_idx = start_idx + 1300-1;
            
        else
            end_idx = triggers(n+1,1)-1;
            
        end
        
    else
        end_idx = start_idx + n_data_points_per_epoch-1;
        
    end
    
    cur_ep_length = end_idx - start_idx +1;
    
    start_end_idx(cur_per,cur_rep,:) = [start_idx end_idx];
    epoched_data(1:cur_ep_length,cur_per,cur_rep,:) = permute(raw_data(:,start_idx:end_idx),[2 3 4 1]);
    
    
    
end


% reshape data

epoch.data = epoched_data;
epoch.start_end = start_end_idx;
epoch.idx = epoch_id;


end

















