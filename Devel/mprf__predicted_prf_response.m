function pred_resp = mprf__predicted_prf_response(model, stimulus, prf,roi)

%%% As it is now, this code is not very flexible, possibly update it, i.e.
%%% change it when necessary. Think about how to implement scrambled
%%% versions or versions with a ranged pRF parameter -> 

if isfield(roi,'idx') && model.params.roi_specific
    
    rois = fieldnames(roi.idx);
    pred_resp = zeros(size(stimulus.im,2), size(roi.mask,1),numel(rois));
    
    for nn = 1:length(rois)
        
        cur_in = roi.idx_out == roi.idx.(rois{nn}).idx;
        
        orig_idx = zeros(size(cur_in));
        orig_idx(cur_in) = find(cur_in);
        
        
        %fprintf('Making predictions for %d surface vertices\n',sum(good_data));
        for n = 1:1000:length(prf.sigma.val)
            
            end_idx = n+999;
            if end_idx > length(prf.sigma.val)
                end_idx = length(prf.sigma.val);
            end
            
            
            cur_good_data = cur_in(n:end_idx);
            cur_orig_idx = orig_idx(n:end_idx);
            
            if prf.sigma.type.fixed
                cur_sigma = prf.sigma.val;
                
            else
                cur_sigma = prf.sigma.val(n:end_idx);
                cur_sigma = cur_sigma(cur_good_data);
            end
            
            if prf.x0.type.fixed
                cur_X0 = prf.x0.val;
                
            else
                cur_X0 = prf.x0.val(n:end_idx);
                cur_X0 = cur_X0(cur_good_data);
                
            end
            
            
            if prf.y0.type.fixed
                cur_Y0 = prf.y0.val;
                
            else
                cur_Y0 = prf.y0.val(n:end_idx);
                cur_Y0 = cur_Y0(cur_good_data);
            end
            
            if prf.beta.type.fixed
                cur_beta = prf.beta.val;
                
            else
                cur_beta = prf.beta.val(n:end_idx);
                cur_beta = cur_beta(cur_good_data);
            end
            
            
            cur_orig_idx = cur_orig_idx(cur_good_data);
            
            RF = rfGaussian2d(stimulus.X,stimulus.Y,cur_sigma,cur_sigma,false, cur_X0,cur_Y0);
            cur_pred = bsxfun(@times, stimulus.im' * RF,cur_beta');
            
            pred_resp(:,cur_orig_idx,nn) = cur_pred;
            
            
        end
        
    end
    
end


end






