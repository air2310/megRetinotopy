function [prf, bs, roi] = mprf__model_get_prf_params(model)
% Load the pRF parameters necessary for the current modeling
% Add roi option here as well, i.e. load ROI mask, combine ROIs in the way
% specified from the GUI.

main_dir = mprf__get_directory('main_dir');
bs_surf_dir = mprf__get_directory('bs_surf_prf_data');
tmp_bs_file = mprf__get_file('bs_headmodel');


if length(tmp_bs_file.fname) > 1
    answer = listdlg('ListString',{tmp_bs_file.fname.name},...
        'PromptString','Multiple head model files found, please select one',...
        'SelectionMode','Single');
else
    answer = 1;
    
end


bs_file = fullfile(tmp_bs_file.fpath, tmp_bs_file.fname(answer).name);

fname_parts = strsplit(tmp_bs_file.fname(answer).name,'_');

tmp = {'pial','white','mid'};

file_prefix = tmp{ceil((find(~cellfun(@isempty, ...
    [strfind(fname_parts,tmp{1}),strfind(fname_parts,tmp{2}),...
    strfind(fname_parts,tmp{3})])))./length(fname_parts))};

file_prefix = [file_prefix '.'];

x_type = get_type(model.params.x0);
sigma_type = get_type(model.params.sigma);
y_type = get_type(model.params.y0);
b_type = get_type(model.params.beta);
ve_type = get_type('varexplained');

beta_fname = get_file_name(b_type);
x_fname = get_file_name(x_type);
y_fname = get_file_name(y_type);
sigma_fname = get_file_name(sigma_type);
ve_fname = get_file_name(ve_type);

beta_fpath = fullfile(main_dir, bs_surf_dir, [file_prefix beta_fname]);
x_fpath = fullfile(main_dir, bs_surf_dir, [file_prefix x_fname]);
y_fpath = fullfile(main_dir, bs_surf_dir, [file_prefix y_fname]);
sigma_fpath = fullfile(main_dir, bs_surf_dir, [file_prefix sigma_fname]);
ve_fpath = fullfile(main_dir, bs_surf_dir, [file_prefix ve_fname]);

prf.x0.val = get_prf_val(x_type, model, x_fpath, 'x0');
prf.y0.val = get_prf_val(y_type, model, y_fpath, 'y0');
prf.sigma.val = get_prf_val(sigma_type, model, sigma_fpath, 'sigma');
prf.beta.val = get_prf_val(b_type, model, beta_fpath, 'beta');
prf.ve.val = get_prf_val(ve_type, model, ve_fpath,'ve');

prf.x0.file = x_fpath;
prf.x0.type = x_type;
prf.y0.file = y_fpath;
prf.y0.type = y_type;
prf.sigma.file = sigma_fpath;
prf.sigma.type = sigma_type;
prf.beta.file = beta_fpath;
prf.beta.type = b_type;
prf.ve.file = ve_fpath;
prf.ve.type = ve_type;

bs.model_file = bs_file;

if model.params.roi_specific || model.params.roi_mask
    
    roi_dir = mprf__get_directory('bs_surf_roi');
    all_roi_fname = 'all_rois';
    roi_mask_fname = 'all_rois_mask';
    
    all_roi_fpath = fullfile(main_dir, roi_dir,[file_prefix all_roi_fname ]);
    roi_mask_fpath = fullfile(main_dir, roi_dir,[file_prefix roi_mask_fname ]);
    
    if model.params.roi_specific
        load(fullfile(main_dir, mprf__get_file('bs_tag_to_idx')));
    end
    
    has_all_roi = any(exist(all_roi_fpath, 'file'));
    has_all_roi_mask = any(exist(roi_mask_fpath,'file'));
    
    % We have both maks and ROI idx
    if has_all_roi && has_all_roi_mask
        
        % We onlye have ROI idx
    elseif has_all_roi && ~has_all_roi_mask
        
        % We only have ROI mask, but want roi specific output, so need to
        % specify ROI idx file
    elseif ~has_all_roi && model.params.roi_specific
        
        [fname, fpath] = uigetfile('*', 'Could not locate all ROI file needed for roi specific output, please locate the file');
        
        if isempty(fname)
            error('Invalid file returned')
        end
        
        uni_idx = unique(read_curv(fullfile(fpath, fname)));
        
        if any(isnan(uni_idx))
            uni_idx = uni_idx(~isnan(uni_idx));
        end
        
        nfields = numel(fieldnames(tag_to_idx));
        
        if numel(uni_idx) < nfields
            warning('Only have indices for %d of the %d roi tags',numel(uni_idx), nfields)
        end
        
        all_roi_fpath = fullfile(fpath, fname);
        has_all_roi = true;
        % We do not have any ROI file, but want to mask
    elseif ~has_all_roi && ~has_all_roi_mask && model.params.roi_mask
        
        
        % Check the amount of unique values in file that is returned to
        % determin the type, mask or index file
        [fname, fpath] = uigetfile('*', 'Could not locate all ROI file or all ROI mask, please locate the file');
        
        if isempty(fname)
            error('Invalid file returned')
        end
        
        uni_idx = unique(read_curv(fullfile(fpath, fname)));
        
        if any(isnan(uni_idx))
            uni_idx = [uni_idx(~isnan(uni_idx)) nan];
        end
        
        if numel(uni_idx) == 2;
            
            
        elseif numel(uni_idx) > 2
            
            
        else
            error('Invalid ROI mask file')
        end
        
        has_all_roi_mask = true;
        roi_mask_fpath = fullfile(fpath, fname);
        
    end
    
    if ~has_all_roi
        all_roi_fpath = '';
        
    end
    
    if ~has_all_roi_mask
        roi_mask_fpath = '';
        
    end
    
    roi_info.has_all_roi = has_all_roi;
    roi_info.has_all_roi_mask = has_all_roi_mask;
    roi_info.all_roi_fpath = all_roi_fpath;
    roi_info.roi_mask_fpath = roi_mask_fpath;
    
    if model.params.roi_specific
        roi_info.tag_to_idx = tag_to_idx;
    end
    
    roi = get_roi_data(model, roi_info);
    
else
    roi.mask = ones(size(prf.x0.val));
    roi.idx_out = ones(size(prf.x0.val));
    
end

end


function type = get_type(param)

type.smoothed = any(strfind(param, 'smoothed'));
type.scrambled = any(strfind(param, 'scrambled'));
type.range = any(strfind(param, 'range'));
type.fixed  = any(strfind(param, 'fixed'));
type.absolute  = any(strfind(param, 'absolute'));
type.proportional  = any(strfind(param, 'proportion'));

if any(strfind(param,'mresp'))
    type.ftype = 'mresp';
    
elseif any(strfind(param,'varexplained'))
    type.ftype = 'varexplained';
    
elseif any(strfind(param,'x')) && ~any(strfind(param,'sigma')) && ...
        ~any(strfind(param,'y')) && ~any(strfind(param,'beta'))
    type.ftype = 'x';
    
elseif any(strfind(param,'y'))
    type.ftype = 'y';
    
elseif any(strfind(param,'beta')) && ~any(strfind(param,'recomp_beta'))
    type.ftype = 'beta';
    
elseif any(strfind(param,'beta')) && any(strfind(param,'recomp_beta'))
    type.ftype = 'recomp_beta';
    
elseif any(strfind(param,'sigma'))
    type.ftype = 'sigma';
    
    
    
else
    error('Not implemented option')
    
    
end





end


function fname = get_file_name(type)

if type.smoothed && ~type.absolute
    fname = [type.ftype '_smoothed'];
    
elseif ~type.smoothed && ~type.absolute
    fname = type.ftype;
    
elseif type.absolute
    fname = '';
    
else
    error('Not implemented option')
    
    
end


end

function val = get_prf_val(type, model, fpath, cur_par)


if type.fixed && type.absolute
    if strcmpi(cur_par, 'beta')
        if model.params.beta_equal_beta
            val = model.params.([cur_par '_fix']);
            
        elseif model.params.beta_equal_pred
            val = 1;
            
        else
            error('Not implemented option')
            
            
        end
        
    else
        val = model.params.([cur_par '_fix']);
        
    end
    
    
elseif type.range && type.absolute
    val = model.params.([cur_par '_range']);
    
else
    if exist(fpath,'file')
        val = read_curv(fpath);
        if type.fixed && type.proportional
            val = val .* model.params.([cur_par '_fix']);
            
            
        elseif type.range && type.proportional
            val = val * model.params.([cur_par '_range']);
            
        elseif type.smoothed || type.scrambled || ...
                strcmpi(type.ftype,'recomp_beta') || ...
                strcmpi(type.ftype,'varexplained')
            
        else
            error('Not implemented option')
            
            
        end
        
    else
        error('Could not find file: %s', fpath)
        
    end
end



end

function roi = get_roi_data(model,roi_info)

has_mask = false;

if model.params.roi_mask
    
    if roi_info.has_all_roi_mask
        tmp_mask = read_curv(roi_info.roi_mask_fpath);
        
        if any(isnan(tmp_mask))
            % Assume nans are outside rois
            tmp_mask = double(~isnan(tmp_mask));
            has_mask = true;
            
        elseif numel(unique(tmp_mask)) == 2;
            % Typically, mask is created by isnan(all_rois), so true for nans,
            % false for non-nans. Assuming this case here...
            tmp_mask = double(~tmp_mask);
            has_mask = true;
            
        elseif numel(unique(tmp_mask)) > 2
            % If there are zeros, assuming they mark non-roi vertices
            if any(tmp_mask == 0)
                tmp_mask = tmp_mask > 0;
                has_mask = true;
                
                
            else
                if roi_info.has_all_roi
                    % Try with roi_idx file instead
                    
                else
                    error('Could not reliably process ROI mask')
                end
            end
        else
            if roi_info.has_all_roi
                % Try with roi_idx file instead
                
            else
                error('Could not reliably process ROI mask')
            end
            
        end
        
        
    elseif roi_info.has_all_roi && ~has_mask
        tmp_mask = read_curv(roi_info.all_roi_fpath);
        
        
        if any(isnan(tmp_mask))
            % Assume nans are outside rois
            tmp_mask = double(~isnan(tmp_mask));
            
            
        elseif numel(unique(tmp_mask)) >= 2
            % If there are zeros, assuming they mark non-roi vertices
            if any(tmp_mask == 0)
                tmp_mask = tmp_mask > 0;
                
            else
                error('Could not reliably process ROI mask')
                
            end
            
            
        else
            error('Could not reliably process ROI mask')
            
            
        end
        
    else
        error('No file found to use for ROI mask')
        
        
    end
    
end
roi.mask = tmp_mask;

if model.params.roi_specific
    if roi_info.has_all_roi
        tmp_roi = read_curv(roi_info.all_roi_fpath);
        roi.idx_out = zeros(size(tmp_roi));
        
        roi_tag = roi_info.tag_to_idx;
        roi_fields = fieldnames(roi_tag);
        
        skip_rois = {};
        
        cnt = 0;
        for n = 1:length(roi_fields)
            cur_roi = roi_fields{n};
            [cur_lr, cur_name, cur_dv] = parse_roi_name(cur_roi);
            
            
            %%% CHECK IF ROI IS IN SKIP_ROI VARIABLE HERE
            if any(cellfun(@(x) strcmp(x, cur_roi),skip_rois))
                
            else
                if model.params.comb_dv_rois || model.params.comb_lr_rois
                    comb_rois = {};
                    these_idx = [];
                    for nn = 1:length(roi_fields)
                        [tmp_lr, tmp_name, tmp_dv] = parse_roi_name(roi_fields{nn});
                        if strcmpi(tmp_name, cur_name)
                            
                            if model.params.comb_dv_rois && ~model.params.comb_lr_rois
                                if strcmpi(tmp_lr, cur_lr)
                                    comb_rois = [comb_rois roi_fields{nn}];
                                    these_idx = [these_idx roi_info.tag_to_idx.(roi_fields{nn})];
                                    skip_rois = [skip_rois roi_fields{nn}];
                                    
                                end
                                
                                
                            elseif ~model.params.comb_dv_rois && model.params.comb_lr_rois
                                if strcmpi(tmp_dv, cur_dv)
                                    comb_rois = [comb_rois roi_fields{nn}];
                                    these_idx = [these_idx roi_info.tag_to_idx.(roi_fields{nn})];
                                    skip_rois = [skip_rois roi_fields{nn}];
                                    
                                end
                                
                                
                            elseif model.params.comb_dv_rois && model.params.comb_lr_rois
                                comb_rois = [comb_rois roi_fields{nn}];
                                these_idx = [these_idx roi_info.tag_to_idx.(roi_fields{nn})];
                                skip_rois = [skip_rois roi_fields{nn}];
                            end
                        end
                    end
                else

                    comb_rois = roi_fields{n};
                    these_idx = roi_info.tag_to_idx.(roi_fields{n});
                    
                end
                
                if model.params.comb_dv_rois && ~model.params.comb_lr_rois
                    if isempty(cur_lr)
                        
                    else
                        cur_name = [cur_lr '_' cur_name];
                    end
                elseif ~model.params.comb_dv_rois && model.params.comb_lr_rois
                    if isempty(cur_dv)
                        
                    else
                        cur_name = [cur_name '_' cur_dv];
                    end
                elseif model.params.comb_dv_rois && model.params.comb_lr_rois
                    
                    
                else
                    cur_name = cur_roi;
                    
                end
                cnt = cnt + 1;
                
                roi.idx_out(ismember(tmp_roi,these_idx)) = cnt;
                roi.idx.(cur_name).idx = cnt;
                roi.idx.(cur_name).orig_idx = these_idx;
                roi.idx.(cur_name).orig_rois = comb_rois;
                
                
            end
        end
        
        
    else
        error('No file found to use for ROI mask')
        
    end
    

    
end





end

function [cur_lr, cur_name, cur_dv] = parse_roi_name(roi_name)

cur_lr = '';
cur_name = '';
cur_dv = '';

tmp = strsplit(roi_name,'_');
if length(tmp) == 1
    
    % Expecting just the name of the ROI:
    if strcmpi(tmp{1},'left') || strcmpi(tmp{1},'right') || ...
            strcmpi(tmp{1},'d') || strcmpi(tmp{1},'v')
        error('Could not parse ROI tag %s', cur_roi)
    else
        
        
    end
    
    cur_name = tmp{1};
    %%% NO NEED TO CONTINUE HERE.... THERE IS NO LEFT OR RIGHT, OR DORSAL OR VENTRAL VERSION OF THIS ROI
    
elseif length(tmp) == 2
    % Expecting either lr_name, or name_dv
    if strcmpi(tmp{1},'left') || strcmpi(tmp{1},'right')
        cur_lr = tmp{1};
        
        if ~strcmpi(tmp{2},'d') && ~strcmpi(tmp{2},'v')
            cur_name = tmp{2};
            
        else
            error('Could not parse ROI tag %s', cur_roi)
        end
        
    elseif ~strcmpi(tmp{1},'left') && ~strcmpi(tmp{1},'right') && ...
            ~strcmpi(tmp{1},'d') && ~strcmpi(tmp{1},'v')
        cur_name = tmp{1};
        
        if strcmpi(tmp{2},'d') || strcmpi(tmp{2},'v')
            cur_dv = tmp{2};
            
            
        else
            error('Could not parse ROI tag %s', cur_roi)
            
        end
        
    else
        error('Could not parse ROI tag %s', cur_roi)
        
    end
    
elseif length(tmp) == 3
    
    if strcmpi(tmp{1},'left') || strcmpi(tmp{1},'right')
        cur_lr = tmp{1};
        
        if strcmpi(tmp{2},'left') || strcmpi(tmp{2},'right') || ...
                strcmpi(tmp{2},'d') || strcmpi(tmp{2},'v')
            error('Could not parse ROI tag %s', cur_roi)
            
        else
            cur_name = tmp{2};
            
        end
        
        if strcmpi(tmp{3},'d') || strcmpi(tmp{3},'v')
            cur_dv = tmp{3};
            
        else
            
            error('Could not parse ROI tag %s', cur_roi)
            
        end
        
    else
        error('Could not parse ROI tag %s', cur_roi)
        
    end
    
else
    error('Could not parse ROI tag %s', cur_roi)
end




end


