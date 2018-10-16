function epoch = mprfPreprocessWrapper(epoch)

varThreshold        = [0.05 20];
badChannelThreshold = 0.2;
badEpochThreshold   = 0.2;
verbose             = true;

first_pos_idx = abs([0; diff(epoch.idx)]) > 0 & abs([0; diff(epoch.idx)]) < 10;
blink_idx = epoch.idx == 20;

remove_epochs = first_pos_idx | blink_idx;

data_reshaped = epoch.data;
orig_dim = size(epoch.data);

epoch = rmfield(epoch,'data');

fprintf('Removing blink and first bar position epochs...\n')
data_reshaped = data_reshaped(:,~remove_epochs,:,:);
data_dim = size(data_reshaped);

fprintf('Done.\n')

fprintf('Reshaping output data...\n')
data_reshaped = reshape(data_reshaped,data_dim(1), [], data_dim(4)); % Assuming reshape does what I expect it to do....
fprintf('Done.\n')

fprintf('Preprocessing data...\n')
[data_reshaped, badChannels, badEpochs]  = nppPreprocessData(data_reshaped, ...
    varThreshold, badChannelThreshold, badEpochThreshold,verbose);
fprintf('Done.\n')

fprintf('Removing bad epochs and channels...\n')
data_reshaped(:,badEpochs,:) = nan;
data_reshaped(:,:,badChannels) = nan;
fprintf('Done.\n')

fprintf('Preallocating output variable.\n Estimated %d bites of memory.\n',...
    prod([orig_dim 8]))

epoch.data = nan(orig_dim);
epoch.data(:,~remove_epochs,:,:) = reshape(data_reshaped,data_dim);
fprintf('Done.\n')

clear data_reshaped

epoch.preproc.badChannels = badChannels;
epoch.preproc.badEpochs = badEpochs;
epoch.preproc.badChannel_thr = badChannelThreshold;
epoch.preproc.badEpoch_thr = badEpochThreshold;
epoch.preproc.var_thr = varThreshold;
epoch.preproc.first_pos_idx = first_pos_idx;
epoch.preproc.blink_idx = blink_idx;




end

