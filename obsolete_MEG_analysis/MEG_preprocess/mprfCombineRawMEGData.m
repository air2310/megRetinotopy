%% mprfCombineRawMEGData

% Directory where the raw data is located:
raw_dir = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wlsubj030/raw/';

% Name of data file to load:
fname = 'R0942_Ret*';

% Load and concatenate time series
data = meg_load_sqd_data(raw_dir, fname);


% Rewrite data as one sqd file where data are NumSamples x NumChannels
source = fullfile(raw_dir, 'R0942_Ret1.sqd');
destination = fullfile(raw_dir, 'R0942_Ret_combined.sqd');

[info, err] = sqdwrite(source, destination, data');