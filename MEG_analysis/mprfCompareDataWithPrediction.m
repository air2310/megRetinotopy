

measure = 'amp_ratio';
add_stimulus = true;
pred_stim = [];
rootpath = mprfRootPath;

cur_dir = pwd;
cd(rootpath(1:strfind(rootpath,'/Code/mprfSession')))
[data_file, data_path] = uigetfile('*','Please select data file to use');
[pred_file, pred_path] = uigetfile('*','Please select prediction to use');
[gain_file, gain_path] = uigetfile('*','Please select predicted gain to use');

% sub_dir = pred_path(1:strfind(pred_path, '/prediction/synthetic/Preprocessing'));
% 
% cd(sub_dir)

fprintf('Loading data from %s\n',data_file)
tmp = load(fullfile(data_path, data_file));
fprintf('Done.\n')

fnames = fieldnames(tmp);
data = tmp.(fnames{1});

fprintf('Loading data from %s\n',pred_file)
tmp = load(fullfile(pred_path, pred_file));
fprintf('Done.\n')

fnames = fieldnames(tmp);
pred = tmp.(fnames{1});

load(fullfile(gain_path, gain_file));

data_pp = mprfExtractAmplitudeMeasureFromTimeSeries(data, measure);
pred_pp = mprfExtractAmplitudeMeasureFromTimeSeries(pred, measure);

data_pp.data.preproc = data.preproc;
data_pp.data.idx = data.idx;
data_pp.start_end = data.start_end;

pred_pp.data.preproc = pred.preproc;
pred_pp.data.idx = pred.idx;
pred_pp.start_end = pred.start_end;
save('Plot_data','data_pp','pred_pp','gain');

%%
periods.blank = [3:5 30:32 57:59 84:86 111:113 138:140];
periods.blink = [1 2 28 29 55 56 82 83 109 110 136 137];
periods.stim = setdiff(1:140,[periods.blink periods.blank]);

if add_stimulus
   [stim_file, stim_path] = uigetfile('*','Please select stimulus used to generate the predictions/acquire the MEG data');
   load(fullfile(stim_path, stim_file))
    
end

mprfPlotMeasuredAndPredictedMEGData(data_pp, pred_pp, periods, pred_stim, measure);

mprfMakeDataPredAmplitudeContrasts(data, pred, pred_stim.full_im, periods);


snr_tmp.r = nan(1,size(data_pp.snr_tmp,2));
snr_data.r = nan(1,size(data_pp.snr_data,2));
snr_median.r = nan(1,size(data_pp.snr_median,2));
snr_median_gain.r = nan(1,size(data_pp.snr_median,2));
snr_co_gain.r = nan(1,size(data_pp.snr_median,2));

snr_tmp.p = nan(1,size(data_pp.snr_tmp,2));
snr_data.p = nan(1,size(data_pp.snr_data,2));
snr_median.p = nan(1,size(data_pp.snr_median,2));
snr_median_gain.p = nan(1,size(data_pp.snr_median,2));
snr_co_gain.p = nan(1,size(data_pp.snr_median,2));

for n = 1:size(data_pp.snr_median,2)
    
    cur_data = data_pp.snr_median(:,n);
    cur_pred = pred_pp.snr_median(:,n);
    
    sel = ~isnan(cur_data) & ~isnan(cur_pred);
    
    [tmp, tmp2] = corrcoef(cur_data(sel), cur_pred(sel));
    
    snr_median.r(n) = tmp(1,2);
    snr_median.p(n) = tmp2(1,2);
    
    cur_data = data_pp.snr_tmp(:,n);
    cur_pred = pred_pp.snr_tmp(:,n);
    
    sel = ~isnan(cur_data) & ~isnan(cur_pred);
    
    [tmp, tmp2] = corrcoef(cur_data(sel), cur_pred(sel));
    
    snr_tmp.r(n) = tmp(1,2);
    snr_tmp.p(n) = tmp2(1,2);
    
    
    
    cur_data = data_pp.snr_data(:,n);
    cur_pred = pred_pp.snr_data(:,n);
    
    sel = ~isnan(cur_data) & ~isnan(cur_pred);
    
    [tmp, tmp2] = corrcoef(cur_data(sel), cur_pred(sel));
    
    snr_data.r(n) = tmp(1,2);
    snr_data.p(n) = tmp2(1,2);
    
    
    cur_data = data_pp.snr_median(:,n);
    cur_pred = abs(gain.M_02(:,n));
    
    sel = ~isnan(cur_data) & ~isnan(cur_pred);
    
    [tmp, tmp2] = corrcoef(cur_data(sel), cur_pred(sel));
    
    snr_median_gain.r(n) = tmp(1,2);
    snr_median_gain.p(n) = tmp2(1,2);
    
    
    cur_data = data_pp.sl_co(:,n);
    cur_pred = abs(gain.M_02(:,n));
    
    sel = ~isnan(cur_data) & ~isnan(cur_pred);
    
    [tmp, tmp2] = corrcoef(cur_data(sel), cur_pred(sel));
    
    snr_co_gain.r(n) = tmp(1,2);
    snr_co_gain.p(n) = tmp2(1,2);
    
    
    cur_data = pred_pp.sl_co(:,n);
    cur_pred = abs(gain.M_02(:,n));
    
    sel = ~isnan(cur_data) & ~isnan(cur_pred);
    
    [tmp, tmp2] = corrcoef(cur_data(sel), cur_pred(sel));
    
    pred_gain_r(n) = tmp(1,2);
    pred_gain_p(n) = tmp2(2,2);
    
    
    
end

fh = figure;
pred_gain_r(isnan(pred_gain_r)) = 0;
fh = megPlotMap(pred_gain_r,[-1 1],fh, jet(256),'Correlation Prediction - Gain');


[snr_median.rs, snr_median.idx] = sort(snr_median.r,'descend');

snr_median.rs_02 = snr_median.rs(~isnan(snr_median.rs));
snr_median.idx_02 = snr_median.idx(~isnan(snr_median.rs));

snr_median.idx_02(1:10)
mprfPlotHeadLayout(snr_median.idx_02(1:10));
title('Median')



snr_median.rs_map = snr_median.r;
snr_median.rs_map(isnan(snr_median.rs_map)) = 0;

fh = figure; hold on;
fh = megPlotMap(snr_median.rs_map,[-1 1],fh, jet(256),'Correlation data - pred, median');



[snr_median_gain.rs, snr_median_gain.idx] = sort(snr_median_gain.r,'descend');

snr_median_gain.rs_02 = snr_median_gain.rs(~isnan(snr_median_gain.rs));
snr_median_gain.idx_02 = snr_median_gain.idx(~isnan(snr_median_gain.rs));

snr_median_gain.idx_02(1:10)
mprfPlotHeadLayout(snr_median_gain.idx_02(1:10));
title('Median gain')



snr_median_gain.rs_map = snr_median_gain.r;
snr_median_gain.rs_map(isnan(snr_median_gain.rs_map)) = 0;

fh = figure; hold on;
fh = megPlotMap(snr_median_gain.rs_map,[-1 1],fh, jet(256),'Correlation data - pred, median');



[snr_co_gain.rs, snr_co_gain.idx] = sort(snr_co_gain.r,'descend');

snr_co_gain.rs_02 = snr_co_gain.rs(~isnan(snr_co_gain.rs));
snr_co_gain.idx_02 = snr_co_gain.idx(~isnan(snr_co_gain.rs));

snr_co_gain.idx_02(1:10)
mprfPlotHeadLayout(snr_co_gain.idx_02(1:10));
title('Coherence gain')



snr_co_gain.rs_map = snr_co_gain.r;
snr_co_gain.rs_map(isnan(snr_co_gain.rs_map)) = 0;

fh = figure; hold on;
fh = megPlotMap(snr_co_gain.rs_map,[-1 1],fh, jet(256),'Correlation data - pred, median');


% 
% [snr_data.rs, snr_data.idx] = sort(snr_data.r,'descend');
% 
% snr_data.rs_02 = snr_data.rs(~isnan(snr_data.rs));
% snr_data.idx_02 = snr_data.idx(~isnan(snr_data.rs));
% 
% snr_data.idx_02(1:10)
% mprfPlotHeadLayout(snr_data.idx_02(1:10));
% title('Data')
% 
% 
% snr_data.rs_map = snr_data.r;
% snr_data.rs_map(isnan(snr_data.rs_map)) = 0;
% 
% fh = figure; hold on;
% fh = megPlotMap(snr_data.rs_map,[-1 1],fh, jet(256),'Correlation data - pred, data');
% 
% 
% 
% [snr_tmp.rs, snr_tmp.idx] = sort(snr_tmp.r,'descend');
% 
% snr_tmp.rs_02 = snr_tmp.rs(~isnan(snr_tmp.rs));
% snr_tmp.idx_02 = snr_tmp.idx(~isnan(snr_tmp.rs));
% 
% snr_tmp.idx_02(1:10)
% mprfPlotHeadLayout(snr_tmp.idx_02(1:10));
% title('Tmp')
% 
% 
% snr_tmp.rs_map = snr_tmp.r;
% snr_tmp.rs_map(isnan(snr_tmp.rs_map)) = 0;
% 
% fh = figure; hold on;
% fh = megPlotMap(snr_tmp.rs_map,[-1 1],fh, jet(256),'Correlation data - pred, tmp');

plot_chnl = [17    26    34    36    42    66    81    97    98   103   104];
%snr_co_gain.idx_02(1:5)
mprfPlotMeasuredAndPredTseries(data_pp.snr_median, pred_pp.snr_median, periods,plot_chnl ,gain)



mean_co = nanmean(data_pp.sl_co);
mean_co(isnan(mean_co)) = 0;

fh = figure;
fh = megPlotMap(mean_co, [0 .25],fh);




