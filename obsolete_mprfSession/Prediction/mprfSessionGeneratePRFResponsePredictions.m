function [pred_resp, pred_stim, stim_file_to_use] = mprfSessionGeneratePRFResponsePredictions(surf, pred_stim,modality, thr)

if ~exist('cellfind','file')
    tbUse('vistasoft')
end

global mprfSESSION



if ~exist('pred_stim','var') || isempty(pred_stim)
    
    stim_file_to_use = mprfSessionGetPredictionStimulus;
    tmp = load(stim_file_to_use);
    fnames = fieldnames(tmp);
    pred_stim = tmp.(fnames{1});
    clear tmp
    
end

if strcmpi(surf.type,'brainstorm')
    if isfield(surf,'name') && ~isempty(surf.name)
        
        
        
    else
        load(fullfile(pwd,'mprfSESSION.mat'))
        [~, surf.name] = fileparts(bs_model.SurfaceFile);
        
    end
    
    fname_parts = strsplit(surf.name,'_');
    
    tmp = {'pial','white','mid'};
    
    file_prefix = tmp{ceil(find(~cellfun(@isempty, [cellfun(@(x) strfind(x,tmp{1}), fname_parts, 'UniformOutput',false) ...
        cellfun(@(x) strfind(x,tmp{2}), fname_parts, 'UniformOutput',false) ...
        cellfun(@(x) strfind(x,tmp{3}), fname_parts, 'UniformOutput',false)])) ...
        ./length(fname_parts))};
    
    file_prefix = [file_prefix '.'];
    
    
    
elseif strcmpi(surf_type,'freesurfer')
    
    
    
end


needed_par = {'sigma_smoothed','x_smoothed','y_smoothed','recomp_beta','varexplained'};

in_dir = dir(fullfile(mprfSESSION.prf_exp.bs_surface_data,[file_prefix '*']));
in_dir_fnames = {in_dir.name};

for n = 1:length(needed_par)
        
    %cellfun(@(x)regexp(x, [file_prefix needed_par{n} '\>']), in_dir_fnames, 'UniformOutput',false)
    
    is_file = ~cellfun(@isempty, strfind(in_dir_fnames, [file_prefix needed_par{n}]));
    
    if sum(is_file) == 1
        cur_path = fullfile(mprfSESSION.prf_exp.bs_surface_data, in_dir_fnames{is_file});
    else
        poss_files = in_dir_fnames(is_file);
        
        answer = listdlg('ListString',poss_files,...
            'PromptString',sprintf('Multiple files for %s have been found, please select one',needed_par{n}),...
            'SelectionMode','single');
        
        tmp = find(is_file);
        tmp = tmp(answer);
        
        cur_path = fullfile(mprfSESSION.prf_exp.bs_surface_data, in_dir(tmp).name);
        
        
    end
    cur_data = read_curv(cur_path);
    
    if strcmpi(needed_par{n},'x_smoothed') || strcmpi(needed_par{n},'x')
        X0 = cur_data;
        
        
    elseif strcmpi(needed_par{n},'y_smoothed') || strcmpi(needed_par{n},'y')
        Y0 = cur_data;
        
        
    elseif strcmpi(needed_par{n},'sigma_smoothed') || strcmpi(needed_par{n},'sigma')
        sigma = cur_data;
        
        
    elseif strcmpi(needed_par{n},'recomp_beta') || ...
            strcmpi(needed_par{n},'beta') || ...
            strcmpi(needed_par{n},'mresp') || ...
            strcmpi(needed_par{n},'mresp_smoothed')
        beta = cur_data;
        
        
        
    elseif strcmpi(needed_par{n},'varexplained')
        ve = cur_data;
        
    else
        
        error('Unknown data type');
        
        
    end
    
    
end

if strcmpi(thr, 'threshold')
    ve_thr = 0.1; % Absolute
    beta_thr = 95; % Percentile
    
    good_data = ve > ve_thr & beta < prctile(beta, beta_thr);
    orig_idx = zeros(size(good_data));
    orig_idx(good_data) = find(good_data);
    
    
elseif strcmpi(thr,'roi')
    good_data = read_curv(fullfile(mprfSESSION.roi_exp.brainstorm,[file_prefix 'all_rois_mask']));
    good_data = ~good_data;
    
    orig_idx = zeros(size(good_data));
    orig_idx(good_data) = find(good_data);
    
end


if strcmpi(modality,'fmri')
    pred_resp = nan(size(pred_stim.im,2), numel(good_data));
    
elseif strcmpi(modality,'meg')
    pred_resp = zeros(size(pred_stim.im,2), numel(good_data));
    
end

fprintf('Making predictions for %d surface vertices\n',sum(good_data));
for n = 1:1000:length(sigma)
    
    end_idx = n+999;
    if end_idx > length(sigma)
        end_idx = length(sigma);
    end
    
    cur_sigma = sigma(n:end_idx);
    cur_X0 = X0(n:end_idx);
    cur_Y0 = Y0(n:end_idx);
    cur_beta = beta(n:end_idx);
    cur_good_data = good_data(n:end_idx);
    cur_orig_idx = orig_idx(n:end_idx);
    
    
    cur_orig_idx = cur_orig_idx(cur_good_data);
    cur_sigma = cur_sigma(cur_good_data);
    cur_X0 = cur_X0(cur_good_data);
    cur_Y0 = cur_Y0(cur_good_data);
    cur_beta = cur_beta(cur_good_data);
    RF = rfGaussian2d(pred_stim.X,pred_stim.Y,cur_sigma,cur_sigma,false, cur_X0,cur_Y0);
    cur_pred = bsxfun(@times, pred_stim.im' * RF,cur_beta');
    
    pred_resp(:,cur_orig_idx) = cur_pred;
    
    
end
fprintf('Done.\n')





end

