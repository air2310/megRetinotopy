function mprfSession_run_model_server(params,meg_data_file_path,plot_figure)

% % %% TO DO:
% % - Add options to the GUI that allow to combine ROIs
% dorsal and ventral, left and right
% - Set the broad band range in an
% intelligent way -> look at Eline's project for defaults
%
%   Allow run model to set some additional thresholds, like VE and beta
%   threshold
%
% Make gui to load model files that shows the most important options
% describing the file
% Make gui to visualize predicted responses on mesh, also adjust
% mprfSessionRenderDataOnBrainstormSurface and other functions to be
% compatible with new code set up
% Adjust MEG data preprocessing script, make gui as well with some
% options??
%
% Create code that replicates the first attempt and expand on this when
% necessary


%%
% global mprfSESSION
% if exist(fullfile(pwd,'mprfSESSION.mat'),'file')
%     load('mprfSESSION.mat');
%     
% else
%     error('Could not find mprfSESSION file. Please run from session folder')
% end

if ~exist('read_curv','file')
    tbUse({'vistasoft','retMEG'})
end

main_dir = mprf__get_directory('main_dir');

% Get the model options
%[model, stimulus, syn, channels] = mprfSession_model_gui;
model = params.model;
stimulus = params.stimulus;
syn = params.syn;
channels = params.channels;

% If we want to do the reliability check, run that one for the type of
% metric selected. See comments in the respective files for more
% explanation
if strcmpi(model.type,'reliability check')
    
    if model.params.reliability_scans || model.params.reliability_split_half
        
        if model.params.do_sl && model.params.do_bb
            fprintf('Running reliability for stimulus locked and broad band signal\n')
            data = [];
            data = mprf__do_reliability_analysis(model,'stimulus_locked',data,meg_data_file_path);
            mprf__do_reliability_analysis(model,'broadband',data,'meg_data_file_path');
            
            
        elseif  model.params.do_sl && ~model.params.do_bb
            fprintf('Running reliability for stimiulus locked only\n')
            mprf__do_reliability_analysis(model,'stimulus_locked',[],'meg_data_file_path');
            
        elseif  ~model.params.do_sl && model.params.do_bb
            fprintf('Running reliability for broad band only\n')
            mprf__do_reliability_analysis(model,'broadband',[],'meg_data_file_path');
      
        end
    end
    
    % If scramble, original, fixed pRF size of pRF size range:
elseif strcmpi(model.type,'scramble pRF parameters') || ...
        strcmpi(model.type,'run original model') || ...
        strcmpi(model.type,'fix prf size') || ...
        strcmpi(model.type, 'prf size range') ||...
        strcmpi(model.type, 'position (x,y) range')
    
    % Get the modeling parameters from the model variable:
    [prf, bs, roi] = mprf__model_get_prf_params(model);
    
    % Check for how many variables we want to use a range:
    n_range = sum([prf.sigma.type.range
        prf.x0.type.range
        prf.y0.type.range
        prf.beta.type.range]);
    
    % Check how many variables we want to scramble:
    n_scramble = sum([prf.sigma.type.scrambled
        prf.x0.type.scrambled
        prf.y0.type.scrambled
        prf.beta.type.scrambled]);
    
    % Check how many variables to be scramble are grouped:
    n_grouped = sum([prf.sigma.type.grouped
        prf.x0.type.grouped
        prf.y0.type.grouped
        prf.beta.type.grouped]);
    
    
    iter = [];
    
    if n_range >= 1 && n_scramble >= 1
        error('Both a range manipulation and a scrambling is not implemented')
        
    elseif n_range == 0 && n_scramble == 0
        
        % If only one range variable
    elseif n_range == 1 && n_scramble == 0
        tmp_fields = fieldnames(prf);
        iter.method = 'range';
        for n = 1:length(tmp_fields)
            if prf.(tmp_fields{n}).type.range
                
                iter.var = tmp_fields{n};
                iter.n = size(prf.(tmp_fields{n}).val,2);
                
            end
        end
        
        % If more than one range variable (example position (x and y))
    elseif n_range >= 1 && n_scramble == 0
        tmp_fields = fieldnames(prf);
        iter.method = 'range';
        nc=1;
        for n = 1:length(tmp_fields)
            if prf.(tmp_fields{n}).type.range
                
                iter.var{nc} = tmp_fields{n};
                iter.n = size(prf.(tmp_fields{n}).val,2);
                nc =nc+1;
            end
        end
        
        % if only 1 scramble parameter
    elseif n_scramble >= 1 && n_range == 0;
        if n_grouped >=1
            iter.method = 'scramble';
            w_iter_var = cell(1,1);
            w_iter_var2 = '';
            nc = 1;
            
            if prf.x0.type.scrambled
                %nc = nc+1;
                w_iter_var{nc} = 'grouped';
                w_iter_var2 = [w_iter_var2 '_' w_iter_var{nc}];
            end
            
            if prf.y0.type.scrambled
                %nc = nc+1;
                w_iter_var{nc} = 'grouped';
                w_iter_var2 = [w_iter_var2 '_' w_iter_var{nc}];
                
            end
            
            if prf.sigma.type.scrambled
                %nc = nc+1;
                w_iter_var{nc} = 'grouped';
                w_iter_var2 = [w_iter_var2 '_' w_iter_var{nc}];
                
            end
            
            if prf.beta.type.scrambled
                %nc = nc+1;
                w_iter_var{nc} = 'grouped';
                w_iter_var2 = [w_iter_var2 '_' w_iter_var{nc}];
                
            end
            
        else
            iter.method = 'scramble';
            w_iter_var = cell(1,n_scramble);
            w_iter_var2 = '';
            nc = 0;
            
            if prf.x0.type.scrambled
                nc = nc+1;
                w_iter_var{nc} = 'x0';
                w_iter_var2 = [w_iter_var2 '_' w_iter_var{nc}];
            end
            
            if prf.y0.type.scrambled
                nc = nc+1;
                w_iter_var{nc} = 'y0';
                w_iter_var2 = [w_iter_var2 '_' w_iter_var{nc}];
                
            end
            
            if prf.sigma.type.scrambled
                nc = nc+1;
                w_iter_var{nc} = 'sigma';
                w_iter_var2 = [w_iter_var2 '_' w_iter_var{nc}];
                
            end
            
            if prf.beta.type.scrambled
                nc = nc+1;
                w_iter_var{nc} = 'beta';
                w_iter_var2 = [w_iter_var2 '_' w_iter_var{nc}];
                
            end
            
        end
    else
        error('Could not determine scramble parameter');
        %end
        
    %elseif n_range > 1
        %error('Range not implemented for multiple variables')
        
    %else
        %error('Unknown option')
    end
    bs = mprf__get_lead_field(bs);
end
%bs = mprf__get_lead_field(bs);

% BASICALLY WE SHOULD DO SOMETHING DIFFERENT IF WE'RE SCRAMBLING HERE.
% THESE FUNCTIONS ACTUALLY STORE THE PREDICTIONS, WHICH IS FINE IF
% THERE ARE ONLY 15 OF THEM, NOT WHEN THERE ARE A THOUSAND...

if strcmpi(model.type,'run original model') || ...
        strcmpi(model.type,'prf size range') || ...
        strcmpi(model.type,'fix prf size') ||...
        strcmpi(model.type,'position (x,y) range') 

    
    pred_resp = mprf__predicted_prf_response(model, stimulus, prf, roi, iter);
    meg_resp = mprf__compute_predicted_meg_response(bs, pred_resp, channels);
    
    cur_time = datestr(now);
    cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
    
    save_dir = mprf__get_directory('model_predictions');
    save(fullfile(main_dir, save_dir, ['model_predictions_' cur_time]),...
        'prf','bs','roi','model','stimulus','pred_resp','syn','meg_resp','channels');
    
    pred.pred_resp = pred_resp;
    pred.prf = prf;
    pred.bs = bs;
    pred.roi = roi;
    pred.model = model;
    pred.stimulus = stimulus;
    pred.syn = syn;
    pred.meg_resp = meg_resp;
    pred.channels = channels;
    pred.cur_time = cur_time;
    
    
    if strcmpi(model.type,'run original model') || ...
            strcmpi(model.type,'fix prf size')
        
        
        
        mprfSession_run_original_model_server(pred,meg_data_file_path,plot_figure);
        
    elseif strcmpi(model.type,'prf size range') 
        
        mprfSession_run_prf_size_range_model_server(pred,meg_data_file_path,plot_figure);
    
    elseif strcmpi(model.type,'position (x,y) range')
        
       mprfSession_run_prf_position_range_model_server(pred,meg_data_file_path,plot_figure);
        
    end
    
    % if we are scrambling pRF parameters:
elseif strcmpi(model.type, 'scramble prf parameters')
    rng('shuffle');
    model.params.do_sl = true;
    model.params.do_bb = false;
    % Run the original model first to get a baseline measure of the
    % quality of the fit
    % Get predicted pRF responses for main model:
    pred_resp = mprf__predicted_prf_response(model, stimulus, prf, roi);
    
    
    n_cores = model.params.n_cores;
    n_it = model.params.n_iterations;
    
    % Get predicted MEG responses for main model
    meg_resp = mprf__compute_predicted_meg_response(bs, pred_resp, channels);
    
    pred.pred_resp = pred_resp;
    pred.model = model;
    pred.meg_resp = meg_resp;
    pred.channels = channels;
    
    data_in.first_iteration = true;
    % Do an intial fit and return the data
    fprintf('Original model without scrambling \n');
    orig_data = mprfSession_run_original_model_server(pred,meg_data_file_path,plot_figure,data_in);
    
    data_in = orig_data;
    data_in = rmfield(data_in,'cur_corr');
    scramble_corr = nan(size(meg_resp{1},2),n_it);
    
    % Now, create different pRF and MEG response predictions based on
    % the scrambled pRF parameter
    if n_cores > 1
        
        
        if isempty(gcp('nocreate'))
            fprintf('No open pool found\n')
        else
            answer = questdlg('An open matlab pool is found. Do you want to close it or run on a single core',...
                'Open Matlab pool found','Close','Single core','Cancel','Close');
            
            switch lower(answer)
                
                case 'close'
                    delete(gcp);
                    
                case 'single core'
                    pred.model.params.n_cores = 1;
                    n_cores = 1;
                case 'cancel'
                    return
            end
        end
        
        
        % As we are scrambling the pRF parameters, we need to apply the VE
        % and beta thresholds to the parameters before we scramble. We do
        % this here. This is only implemented when assessing all ROIs
        % simultaneously, not for ROI specific models
        
        if ~isfield(roi,'idx') && ~model.params.roi_specific
            cur_in = logical(roi.mask); % These parameters are included
            
        else
            error('Scrambling not implemented for separate ROIs');
            
        end
        % If we need a beta threshold, apply it to the selection. Note that
        % the threshold itself takes into account all the pRF parameters
        % and not just the ones inside the ROI mask (i.e. the selected
        % ones). This is relevant when it comes to the beta threshold, as
        % it is a percentile and not an absolute value. You may want to
        % change this later one. As this code is copied from
        % mprf__predicted_prf_response.m, you may want to make the same
        % changes there as well.
        if model.params.beta_thr
            if model.params.beta_thr_vals(1) == 0
                range = [prctile(prf.beta.val,0) - 1 prctile(prf.beta.val,model.params.beta_thr_vals(2))];
                
            else
                range = [prctile(prf.beta.val,model.params.beta_thr_vals(1)) prctile(prf.beta.val,model.params.beta_thr_vals(2))];
                
            end
            cur_in = cur_in & prf.beta.val > range(1) & prf.beta.val < range(2);
            model.params.beta_thr = false;
        end
        
        if model.params.ve_thr
            range = [model.params.ve_thr_vals(1) model.params.ve_thr_vals(2)];
            cur_in = cur_in & prf.ve.val > range(1) & prf.ve.val < range(2);
            model.params.ve_thr = false;
            
        end
        % Update the roi mask here
        roi.mask = double(cur_in);
        
        % Now use the mask for scrambling:
        roi_idx = find(roi.mask);
        
        configCluster
        c = parcluster;
        mpool = c.parpool(n_cores);
        pctRunOnAll warning('off','all');
        %runParallelCode()
        
        parfor this_it = 1:n_it
            this_prf = prf;
            this_data_in = data_in;
            this_pred = pred;
            
            cur_idx = roi_idx(randperm(size(roi_idx,1)));
            print_str = 'Scrambling:';
            
            for this_var = 1:length(w_iter_var)
                cur_iter_var = w_iter_var{this_var};
                switch lower(cur_iter_var)
                    
                    case 'sigma'
                        this_prf.sigma.val(roi_idx) = this_prf.sigma.val(cur_idx);
                        print_str = [print_str ' Sigma'];
                    case 'y0'
                        this_prf.y0.val(roi_idx) = this_prf.y0.val(cur_idx);
                        print_str = [print_str ' Y0'];
                    case 'x0'
                        this_prf.x0.val(roi_idx) = this_prf.x0.val(cur_idx);
                        print_str = [print_str ' X0'];
                    case 'beta'
                        this_prf.beta.val(roi_idx) = this_prf.beta.val(cur_idx);
                        print_str = [print_str ' beta'];
                    case 'grouped'
                        this_prf.sigma.val(roi_idx) = this_prf.sigma.val(cur_idx);
                        this_prf.y0.val(roi_idx) = this_prf.y0.val(cur_idx);
                        this_prf.x0.val(roi_idx) = this_prf.x0.val(cur_idx);
                        this_prf.beta.val(roi_idx) = this_prf.beta.val(cur_idx);
                        print_str = [print_str ' grouped'];
                    otherwise
                        error('Did not recognize scramble parameter');
                end
                
                
                
            end
            fprintf([print_str '.\n']);
            this_pred_resp = mprf__predicted_prf_response(model, stimulus, this_prf, roi);
            
            this_meg_resp = mprf__compute_predicted_meg_response(bs, this_pred_resp, channels);
            this_pred.meg_resp = this_meg_resp;
            
            
            this_data_in = mprfSession_run_original_model_server(this_pred,meg_data_file_path,plot_figure,this_data_in);
            fprintf('SCRAMBLING ITERATION: %d out of %d \n',this_it,n_it);
            scramble_corr(:,this_it) = this_data_in.cur_corr;
            this_data_in = rmfield(this_data_in,'cur_corr');
            
        end
        if isempty(gcp('nocreate'))
            fprintf('?? no open pool found ??');
        else
            delete(gcp);
        end
        
    elseif n_cores == 1
        
        % As we are scrambling the pRF parameters, we need to apply the VE
        % and beta thresholds to the parameters before we scramble. We do
        % this here. This is only implemented when assessing all ROIs
        % simultaneously, not for ROI specific models
       
        if ~isfield(roi,'idx') && ~model.params.roi_specific
            cur_in = logical(roi.mask); % These parameters are included
            
        else
            error('Scrambling not implemented for separate ROIs');
            
        end
        % If we need a beta threshold, apply it to the selection. Note that
        % the threshold itself takes into account all the pRF parameters
        % and not just the ones inside the ROI mask (i.e. the selected
        % ones). This is relevant when it comes to the beta threshold, as
        % it is a percentile and not an absolute value. You may want to
        % change this later one. As this code is copied from
        % mprf__predicted_prf_response.m, you may want to make the same
        % changes there as well.
        if model.params.beta_thr
            if model.params.beta_thr_vals(1) == 0
                range = [prctile(prf.beta.val,0) - 1 prctile(prf.beta.val,model.params.beta_thr_vals(2))];
                
            else
                range = [prctile(prf.beta.val,model.params.beta_thr_vals(1)) prctile(prf.beta.val,model.params.beta_thr_vals(2))];
                
            end
            cur_in = cur_in & prf.beta.val > range(1) & prf.beta.val < range(2);
            model.params.beta_thr = false;
        end
        
        if model.params.ve_thr
            range = [model.params.ve_thr_vals(1) model.params.ve_thr_vals(2)];
            cur_in = cur_in & prf.ve.val > range(1) & prf.ve.val < range(2);
            model.params.ve_thr = false;
            
        end
        % Update the roi mask here
        roi.mask = double(cur_in);
        
        % Now use the mask for scrambling:
        roi_idx = find(roi.mask);
        
        for this_it = 1:n_it
            
            this_prf = prf;
            this_data_in = data_in;
            this_pred = pred;
            
            cur_idx = roi_idx(randperm(size(roi_idx,1)));
            print_str = 'Scrambling:';
            
            for this_var = 1:length(w_iter_var)
                cur_iter_var = w_iter_var{this_var};
                switch lower(cur_iter_var)
                    
                    case 'sigma'
                        this_prf.sigma.val(roi_idx) = this_prf.sigma.val(cur_idx);
                        print_str = [print_str ' Sigma'];
                    case 'y0'
                        this_prf.y0.val(roi_idx) = this_prf.y0.val(cur_idx);
                        print_str = [print_str ' Y0'];
                    case 'x0'
                        this_prf.x0.val(roi_idx) = this_prf.x0.val(cur_idx);
                        print_str = [print_str ' X0'];
                    case 'beta'
                        this_prf.beta.val(roi_idx) = this_prf.beta.val(cur_idx);
                        print_str = [print_str ' beta'];
                    case 'grouped'
                        this_prf.sigma.val(roi_idx) = this_prf.sigma.val(cur_idx);
                        this_prf.y0.val(roi_idx) = this_prf.y0.val(cur_idx);
                        this_prf.x0.val(roi_idx) = this_prf.x0.val(cur_idx);
                        this_prf.beta.val(roi_idx) = this_prf.beta.val(cur_idx);
                        print_str = [print_str ' grouped'];
                    otherwise
                        error('Did not recognize scramble parameter');
                end

            end
            fprintf([print_str '.\n']);
            this_pred_resp = mprf__predicted_prf_response(model, stimulus, this_prf, roi);
            
            this_meg_resp = mprf__compute_predicted_meg_response(bs, this_pred_resp, channels);
            this_pred.meg_resp = this_meg_resp;
            
            fprintf('SCRAMBLING ITERATION: %d out of %d \n',this_it,n_it);
            this_data_in = mprfSession_run_original_model_server(this_pred,meg_data_file_path,plot_figure,this_data_in);
            
%=======
            %this_data_in = mprfSession_run_original_model_server(this_pred,meg_data_file_path,plot_figure,this_data_in);
%>>>>>>> 847966ef8d7fa9e858d5ad259618662ad511989a
            scramble_corr(:,this_it) = this_data_in.cur_corr;
            this_data_in = rmfield(this_data_in,'cur_corr');
            
        end
    end

    phase_fit_type = [];
    phase_fit_loo = [];
    if strcmp(model.params.metric,'phase ref amplitude') 
        phase_fit_type = model.params.phase_fit;
        if model.params.phase_fit_loo == 1
            phase_fit_loo = 'loo';
        end
    end
    res_dir = mprf__get_directory('model_results');
    main_dir = mprf__get_directory('main_dir');
    rel_dir = 'scramble_prfs';
    cur_time = mprf__get_cur_time;
    
    save_dir = fullfile(main_dir, res_dir, rel_dir, ['Run_stimulus_locked_' w_iter_var2 '_' phase_fit_type '_' phase_fit_loo '_' cur_time]);
    mkdir(save_dir);
    
    results.orig_corr = orig_data.cur_corr;
    results.scrambled_corr = scramble_corr;
    
    save(fullfile(save_dir, 'Results'),'results','model');
    
    
end







if syn.do
    out_name = mprf__make_synthetic_data_set(syn, meg_resp, cur_time,channels);
    
end

if syn.pp
    mprfSession_preprocess_meg_data(true,out_name);
    
end
