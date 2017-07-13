



pred = mprf__load_model_predictions;
[fname, fpath] = uigetfile('*.mat', 'Select raw data to model');
% 
tmp = load(fullfile(fpath, fname));
var_name = fieldnames(tmp);
epoched_data = tmp.(var_name{1});
clear tmp

periods.blank = [3:5 30:32 57:59 84:86 111:113 138:140];
periods.blink = [1 2 28 29 55 56 82 83 109 110 136 137];
periods.stim = setdiff(1:140,[periods.blink periods.blank]);

mprf__coranal_on_meg_data
















