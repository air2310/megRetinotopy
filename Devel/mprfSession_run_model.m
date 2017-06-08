% % %% TO DO:
% % - Add options to the GUI that allow to combine ROIs
% dorsal and ventral, left and right
% - Set the broad band range in an
% intelligent way -> look at Eline's project for defaults
%
%   Allow run model to set some additional thresholds, like VE and beta
%   threshold
%


%%


global mprfSESSION




if exist(fullfile(pwd,'mprfSESSION.mat'),'file')
    load('mprfSESSION.mat');
    
else
    error('Could not find mprfSESSION file. Please run from session folder')
end


[model, stimulus, syn_data] = mprfSession_model_gui;
[prf, bs, roi] = mprf__model_get_prf_params(model);

bs = mprf__get_lead_field(bs);

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


%all_rois = setdiff(unique(roi.mask),0);

% SEL = some_threshold
sel = true(size(roi.mask));

for nn = 1:length(all_rois)
    
    cur_in = roi.mask == all_rois(nn) & sel;
    
    orig_idx = zeros(size(cur_in));
    orig_idx(cur_in) = find(cur_in);
    
    
    %fprintf('Making predictions for %d surface vertices\n',sum(good_data));
    for n = 1:1000:length(prf.sigma.val)
        
        end_idx = n+999;
        if end_idx > length(prf.sigma.val)
            end_idx = length(prf.sigma.val);
        end
        
        cur_sigma = prf.sigma.val(n:end_idx);
        cur_X0 = prf.x0.val(n:end_idx);
        cur_Y0 = prf.y0.val(n:end_idx);
        cur_beta = prf.beta.val(n:end_idx);
        cur_good_data = cur_in(n:end_idx);
        cur_orig_idx = orig_idx(n:end_idx);
        
        
        cur_orig_idx = cur_orig_idx(cur_good_data);
        cur_sigma = cur_sigma(cur_good_data);
        cur_X0 = cur_X0(cur_good_data);
        cur_Y0 = cur_Y0(cur_good_data);
        cur_beta = cur_beta(cur_good_data);
        RF = rfGaussian2d(stimulus.X,stimulus.Y,cur_sigma,cur_sigma,false, cur_X0,cur_Y0);
        cur_pred = bsxfun(@times, stimulus.im' * RF,cur_beta');
        
        pred_resp(:,cur_orig_idx,nn) = cur_pred;
        
        
    end
    
    
    
    
    
end



switch lower(model.type)
    
    
    case 'run original model'
        
        
    case 'equal weight'
        
        
    case 'scramble prf parameters' % NEED SCRAMBLE ITERATIONS
        
        
    case 'fix prf size'
        
        
    case 'pr size range'
        
        
    case 'reliability check'
        
        
        
    case 'fit separate roi predictions'
        
        
        
    otherwise
        
        error('Model type not recognized')
        
        
end






