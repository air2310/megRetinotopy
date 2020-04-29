function mprfSessionComputeGainMatrix

thr_type = 'roi';

load(fullfile(pwd,'mprfSESSION.mat'));
global mprfSESSION

load(mprfSESSION.source.bs_head_model);

[~, surf.name] = fileparts(bs_model.SurfaceFile);
surf.type = 'brainstorm';

[pred_resp, pred_stim, stim_file] = mprfSessionGeneratePRFResponsePredictions(surf, [],'meg',thr_type);
[~, stim_file_name, ~] = fileparts(stim_file);

G = bst_gain_orient(bs_model.Gain,bs_model.GridOrient);
G2 = G(~isnan(G(:,1)),:); % Remove NaNs....
keep_sensors = ~isnan(G(:,1));

M = zeros(size(pred_stim.im,2), size(G,1));
M(:,keep_sensors) = pred_resp * G2';


M_02 = M(:,1:157);

cur_time = datestr(now);
cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';

fname = mprfExportDataPath('pred_gain',['Gain_' cur_time '_' stim_file_name]);

gain.M = M;
gain.M_02 = M_02;
gain.stim = pred_stim;
gain.stim_file = stim_file;

save(fname,'gain');

end


















