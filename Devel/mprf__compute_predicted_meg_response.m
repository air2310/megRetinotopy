function meg_resp = mprf__compute_predicted_meg_response(bs, pred_resp,channels)
% EXPAND THIS FUNCTION TO INCORPORATE DIFFERENT MODEL COMBINATIONS

if size(pred_resp,3) == 1
meg_resp = pred_resp * bs.lead_field(bs.keep_sensors,:)';
meg_resp = meg_resp(:,channels.data+1);

else
    error('Not implemented')
end




end
