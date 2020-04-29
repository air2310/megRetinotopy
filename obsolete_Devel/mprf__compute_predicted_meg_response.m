function meg_resp = mprf__compute_predicted_meg_response(bs, pred_resp,channels)
% EXPAND THIS FUNCTION TO INCORPORATE DIFFERENT MODEL COMBINATIONS

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
