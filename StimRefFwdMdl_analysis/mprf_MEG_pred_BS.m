function [prf,model] = mprf_MEG_pred_BS(subjid)
%% Make predictions for every vertex in brainstorm surface
% Make predictions for every brainstorm vertex
% pred_resp = mprf__predicted_prf_response(model, stimulus, prf, roi, iter);
%
% input - model (information about the different options that describes the model)
%       - MEG stimulus
%       - prf parameters on the brainstorm space
%       - roi vertices on the brainstorm surfaces if available
%
% output - predicted response time series for every vertex on brainstorm
%          surface

% Folder to save output (predicted response at brainstorm vertices to MEG
% stimulus
bs_pred_resp = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Quality_check/%s/pred_resp_bs',subjid);

% MEG stimulus
% CAN ALSO BE OBTAINED FROM MEG DATA FOLDER
templ_file = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Data/MEG/%s/stimFiles/MEG_retinotopy_stimulus_run_1.mat',subjid);
grid_file = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Data/MEG/%s/Stimulus/Other_files/MEG_grid.mat',subjid);

load(templ_file);
load(grid_file);

x_range = [min(grid.Xd(:)) max(grid.Xd(:))];
y_range = [min(grid.Yd(:)) max(grid.Yd(:))];

stim_range = double([min(stimulus.images(:)) max(stimulus.images(:))]);
tmp = stimulus.images(:,:,1);
bk = double(mode(tmp(:)));

new_period = stimulus.seq(stimulus.trigSeq > 0);

im_out_size = [101 101];
meg_stim.resizedIm = zeros([im_out_size, length(new_period)]);

srate = max(grid.Xd(:)) ./ floor(max(im_out_size) / 2);

for n = 1:length(new_period)
    tmp_im = double(stimulus.images(:,:,new_period(n)));
    tmp_im = imresize(tmp_im,im_out_size,'nearest');
    meg_stim.resizedIm(:,:,n) = double(ceil(abs(tmp_im - bk) ./ max(abs(stim_range - bk)))).*(srate.^2);
    
end

[meg_stim.resizedX, meg_stim.resizedY] = meshgrid(linspace(x_range(1), x_range(2), im_out_size(1)),...
    linspace(y_range(1), y_range(2), im_out_size(2)));

keep_window = sum(meg_stim.resizedIm,3) > 0;

meg_stim.im = reshape(meg_stim.resizedIm,[],size(meg_stim.resizedIm,3));
meg_stim.im =  meg_stim.im(keep_window(:),:);
meg_stim.X = meg_stim.resizedX(:);
meg_stim.Y = meg_stim.resizedY(:);
meg_stim.X = meg_stim.X(keep_window);
meg_stim.Y = meg_stim.Y(keep_window);


% stim_file = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Subject_sessions/%s/stimuli/meg/imported_stimulus/meg_stimulus.mat',subjid);
% tmp_stimulus = load(stim_file);
% stimulus  = tmp_stimulus.meg_stim;

% brainstorm surface files folder 
bs_surf_data = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Quality_check/%s/prf_data/surface/brainstorm',subjid);

% describe the model
model_type = 'original (phase ref amplitude) (model fit) (leave one out)';
model = mprf_generate_model_params(model_type);

% path for the prf parameters on the brainstorm surface
fpath_x0 = fullfile(bs_surf_data,'pial.x_smoothed');
fpath_y0 =  fullfile(bs_surf_data,'pial.y_smoothed');
fpath_sigma = fullfile(bs_surf_data,'pial.sigma_smoothed');
fpath_beta = fullfile(bs_surf_data,'pial.recomp_beta');
fpath_ve = fullfile(bs_surf_data,'pial.varexplained');

% pRF parameters on brainstorm surface
prf.x0.val = read_curv(fpath_x0);
prf.y0.val = read_curv(fpath_y0);
prf.sigma.val = read_curv(fpath_sigma);
prf.beta.val = read_curv(fpath_beta);
prf.ve.val = read_curv(fpath_ve);

% predicted pRF response time series on brainstorm surface

% Not using any roi mask. Should be changed to include ROIs (Wang atlas)
roi.mask = ones(size(prf.x0.val));
roi.idx_out = ones(size(prf.x0.val));

rois = 1;
pred_resp = {zeros(size(meg_stim.im,2), size(roi.mask,1),numel(rois))};

nn = 1;% number of rois 
nnn = 1;% number of bootstrapping iterations
cur_in = logical(roi.mask);


if model.params.beta_thr
    if model.params.beta_thr_vals(1) == 0
        range = [prctile(prf.beta.val,0) - 1 prctile(prf.beta.val,model.params.beta_thr_vals(2))];
        
    else
        range = [prctile(prf.beta.val,model.params.beta_thr_vals(1)) prctile(prf.beta.val,model.params.beta_thr_vals(2))];
        
    end
    cur_in = cur_in & prf.beta.val > range(1) & prf.beta.val < range(2);
    
end

if model.params.ve_thr
    range = [model.params.ve_thr_vals(1) model.params.ve_thr_vals(2)];
    cur_in = cur_in & prf.ve.val > range(1) & prf.ve.val < range(2);
end

orig_idx = zeros(size(cur_in));
orig_idx(cur_in) = find(cur_in);
for n = 1:1000:size(prf.sigma.val,1)
    
    end_idx = n+999;
    if end_idx > size(prf.sigma.val,1)
        end_idx = size(prf.sigma.val,1);
    end
    
    
    cur_good_data = cur_in(n:end_idx);
    cur_orig_idx = orig_idx(n:end_idx);
    
    cur_sigma = prf.sigma.val(n:end_idx);
    cur_sigma = cur_sigma(cur_good_data);
    
    
    cur_X0 = prf.x0.val(n:end_idx);
    cur_X0 = cur_X0(cur_good_data);
    
    
    cur_Y0 = prf.y0.val(n:end_idx);
    cur_Y0 = cur_Y0(cur_good_data);
    
    cur_beta = prf.beta.val(n:end_idx);
    cur_beta = cur_beta(cur_good_data);
    
    cur_orig_idx = cur_orig_idx(cur_good_data);
    
    RF = rfGaussian2d(meg_stim.X,meg_stim.Y,cur_sigma,cur_sigma,false, cur_X0,cur_Y0);
    cur_pred = bsxfun(@times, meg_stim.im' * RF,cur_beta');
    
    pred_resp{nnn}(:,cur_orig_idx,nn) = cur_pred;
    
end

disp(sum(sum(isnan(pred_resp{1}))));

% Plot predicted response BS surface
figure, plot(nanmean(pred_resp{1},2))

% save predicted response to MEG stimuli for every brainstorm vertex
save(fullfile(bs_pred_resp,'/pred_resp_bs'),'pred_resp');
