function meg_resp = mprf_MEG_pred_Msen(subjid)

%% Predicted response at MEG sensors to MEG stimuli
% Multiply with gain matrix to get the predictions per vertex
% meg_resp = mprf__compute_predicted_meg_response(bs, pred_resp, channels);
% input - gain matrix
%       - predicted response time series for every vertex on brainstorm
%         surface 
%       - channel information 
%
% output - predicted response time series for every channel on MEG surface

% head model
% CAN ALSO BE OBTAINED FROM BRAINSTORM FOLDER DIRECTLY
bs.model_file = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Data/Brainstorm_db/data/%s/R0774_RETMEG_Block1_5.12.17/headmodel_surf_os_meg.mat',subjid);
bs_model = load(bs.model_file);
bs.lead_field = bst_gain_orient(bs_model.Gain,bs_model.GridOrient);
bs.lead_field2 = bs.lead_field(~isnan(bs.lead_field(:,1)),:); % Remove NaNs....
bs.keep_sensors = ~isnan(bs.lead_field(:,1));

% predicted response at brainstorm vertices to MEG stimulus
bs_pred_resp = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Quality_check/%s/pred_resp_bs',subjid);
load(fullfile(bs_pred_resp,'pred_resp_bs.mat'));

% MEG  Channels
channels.data = 0:156; % MEG channels in which MEG data is acquired
channels.trigger = 160:167; % MEG channels used as triggers
channels.diode = 191; % Channel for collecting diode values



if iscell(pred_resp)
    sz_cell = size(pred_resp,2);
    meg_resp = cell(size(pred_resp));
    
    for n = 1:sz_cell
        fprintf('Iteration: %d\n',n)
        cur_pred = pred_resp{n};
        
        
        if size(cur_pred,3) == 1
            meg_resp{n} = cur_pred * bs.lead_field(bs.keep_sensors,:)';
            meg_resp{n} = meg_resp{n}(:,channels.data+1);
            
        else
            % loop over ROIs here
            error('Not implemented')
        end
        
    end
end   