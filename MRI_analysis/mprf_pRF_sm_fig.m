function mprf_pRF_sm_fig(dirPth)
% plots after mprf_pRF_sm
% 1) distribution of variance explained values for all voxels
% 2) prf size vs eccentricity (different rois)
%                             (original and smoothed)
% 3) surface plots (original and smoothed)

% File paths
% ----------
rootDir = dirPth.rootPth;
prf_dir_mrv = dirPth.fmri.saveDataPth_prfMrv;
prf_data_mrVmat = strcat(rootDir,prf_dir_mrv(2:end),'/mat'); % original and smoothed pRF parameters in mrVista space in .mat
mrSession_dir = dirPth.fmri.mrvPth; 
coords_path = strcat(rootDir,mrSession_dir(2:end),'/Gray/');

anat_dir = dirPth.mri.anatPth; 
roi_path = strcat(rootDir,anat_dir(2:end),'/ROIs/');
% ----------

data_file = dir(fullfile(prf_data_mrVmat,'*.mat'));
load(fullfile(prf_data_mrVmat,data_file.name));

% original prf parameters
ve = prf_par_exp.varexplained;
sigma = prf_par_exp.sigma;
x = prf_par_exp.x; y = prf_par_exp.y;
[pol,ecc] = cart2pol(x,y);
pol= mod(pol,2*pi);
beta = prf_par_exp.beta;

% smoothed prf parameters
sigma_sm = prf_par_exp.sigma_smoothed;
x_sm = prf_par_exp.x_smoothed; y_sm = prf_par_exp.y_smoothed;
[pol_sm,ecc_sm] = cart2pol(x_sm,y_sm);
pol_sm= mod(pol_sm, 2*pi);
recomp_beta = prf_par_exp.recomp_beta;

% histogram of variance explained 
%--------------------------------
figure(1); set(gcf, 'Color', 'w', 'Position', [102   999   1920   400])
hist(ve,100);
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
title('variance explained of prf fits');
xlabel('variance explained');

% prf size vs eccentricity
%-------------------------
% All voxels
figure(2),set(gcf, 'Color', 'w', 'Position', [10   10   1920/3   1080/2])
sm_mask = ve > 0;
ve_thr_mask = ve > 0.4;
c = sm_mask;
scatter(ecc,sigma,[],c);
title('voxels used for smoothing (ve>0)');
xlabel('eccentricity');
ylabel('prf size');

% fits for different visual areas


rois = [{'WangAtlas_V1_Combined'},{'WangAtlas_V2_Combined'},{'WangAtlas_V3_Combined'},...
        {'WangAtlas_IPS_Combined'},{'WangAtlas_LO_Combined'},{'WangAtlas_V3AB_Combined'}...
        {'WangAtlas_MT_Combined'},{'WangAtlas_VO1'},{'WangAtlas_VO2'},{'WangAtlas_FEF'},...
        {'WangAtlas_hV4'}];
    
rois = rois';    
    
num_roi = length(rois);
roi_fname = cell(num_roi,1);

for roi_idx = 1:num_roi
    roi_fname{roi_idx,1} = fullfile(roi_path,strcat(rois{roi_idx},'.mat'));
end    

ROI_params = table(rois,roi_fname);

% Load coordinate file
coordsFile = fullfile(coords_path,'coords.mat');
load(coordsFile);


% Apply these thresholds on the pRF parameters for both the conditions
model_data_thr = cell(num_roi,1);

Var_Exp_Thr = 0.4;
Ecc_Thr = [0 10];
for roi_idx = 1:num_roi
    %Load the current roi
    load(ROI_params.roi_fname{roi_idx});
    
    % find the indices of the voxels from the ROI intersecting with all the voxels
    [~, indices_mean] = intersect(coords', ROI.coords', 'rows' );
    
    model_data = cell(1);
    % Current model parameters- contains x,y, sigma,
    model_data{1}.varexp = ve(:,indices_mean);
    model_data{1}.ecc = ecc(:,indices_mean);
    model_data{1}.sigma = sigma(:,indices_mean);
    model_data{1}.pol = pol(:,indices_mean);
    model_data{1}.x = x(:,indices_mean);
    model_data{1}.y = y(:,indices_mean);
    model_data{1}.beta = beta(:,indices_mean);
    
    model_data{1}.ecc_sm = ecc_sm(:,indices_mean);
    model_data{1}.sigma_sm = sigma_sm(:,indices_mean);
    model_data{1}.pol_sm = pol_sm(:,indices_mean);
    model_data{1}.x_sm = x_sm(:,indices_mean);
    model_data{1}.y_sm = y_sm(:,indices_mean);
    model_data{1}.recomp_beta = recomp_beta(:,indices_mean);
    
    % preallocate variables
    index_thr_tmp = cell(1);
       
    index_thr_tmp{1} = model_data{1}.varexp > Var_Exp_Thr & model_data{1}.ecc < Ecc_Thr(2) & model_data{1}.ecc > Ecc_Thr(1) ;
        
    % Determine the thresholded indices for each of the ROIs
    roi_index{roi_idx,1} = index_thr_tmp{1} ;
       
    % Current model parameters- contains x,y, sigma,
    model_data_thr{roi_idx,1}.varexp = model_data{1}.varexp(roi_index{roi_idx,1});
    model_data_thr{roi_idx,1}.ecc = model_data{1}.ecc(roi_index{roi_idx,1});
    model_data_thr{roi_idx,1}.sigma = model_data{1}.sigma(roi_index{roi_idx,1});
    model_data_thr{roi_idx,1}.pol = model_data{1}.pol(roi_index{roi_idx,1});
    model_data_thr{roi_idx,1}.x = model_data{1}.x(roi_index{roi_idx,1});
    model_data_thr{roi_idx,1}.y = model_data{1}.y(roi_index{roi_idx,1});
    model_data_thr{roi_idx,1}.beta = model_data{1}.beta(roi_index{roi_idx,1});
    
    model_data_thr{roi_idx,1}.ecc_sm = model_data{1}.ecc_sm(roi_index{roi_idx,1});
    model_data_thr{roi_idx,1}.sigma_sm = model_data{1}.sigma_sm(roi_index{roi_idx,1});
    model_data_thr{roi_idx,1}.pol_sm = model_data{1}.pol_sm(roi_index{roi_idx,1});
    model_data_thr{roi_idx,1}.x_sm = model_data{1}.x_sm(roi_index{roi_idx,1});
    model_data_thr{roi_idx,1}.y_sm = model_data{1}.y_sm(roi_index{roi_idx,1});
    model_data_thr{roi_idx,1}.recomp_beta = model_data{1}.recomp_beta(roi_index{roi_idx,1});
    
end

% Store the thresholded pRF values in a table
add_t_1 = table(model_data_thr, 'VariableNames',{'prf_params'});
ROI_params = [ROI_params add_t_1];

figure(4),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
c = [0.5 1 0]; %[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
c_sm = [0 0.5 1]; % ;[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
suptitle('pRF size vs ecc')
for roi_idx = 1:num_roi
    subplot(6,2,roi_idx);
    scatter(ROI_params.prf_params{roi_idx}.ecc,ROI_params.prf_params{roi_idx}.sigma,[],c,'*');  hold on;
    scatter(ROI_params.prf_params{roi_idx}.ecc_sm,ROI_params.prf_params{roi_idx}.sigma_sm,[],c_sm,'*');
    
    tmp = strsplit(ROI_params.rois{roi_idx},'_');
    legend([tmp(2) strcat(tmp(2),'sm')],'Location','NorthWest')
    legend('Location','NorthWest')
    ylim([0 15]);
end

% Sigma histogram
%-----------------------
figure(6),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Sigma')
for roi_idx = 1:num_roi
    subplot(6,2,roi_idx);
    hist(ROI_params.prf_params{roi_idx}.sigma,100)
    tmp = strsplit(ROI_params.rois{roi_idx},'_');
    title(tmp(2))   
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end

% Smoothed sigma histogram
%-----------------------
figure(7),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Sigma smoothed')
for roi_idx = 1:num_roi
    subplot(6,2,roi_idx);
    hist(ROI_params.prf_params{roi_idx}.sigma_sm,100)
    tmp = strsplit(ROI_params.rois{roi_idx},'_');
    title(tmp(2))   
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end


% Beta histogram
%-----------------------
figure(8),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Beta')
for roi_idx = 1:num_roi
    subplot(6,2,roi_idx);
    hist(ROI_params.prf_params{roi_idx}.beta,100)
    tmp = strsplit(ROI_params.rois{roi_idx},'_');
    title(tmp(2))   
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end

% recomputed beta histogram
%-----------------------
figure(9),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Recomputed beta')
for roi_idx = 1:num_roi
    subplot(6,2,roi_idx);
    hist(ROI_params.prf_params{roi_idx}.recomp_beta,100)
    tmp = strsplit(ROI_params.rois{roi_idx},'_');
    title(tmp(2))   
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end

figure(10),set(gcf, 'Color', 'w', 'Position', [10   10   1920/2   1080]);
c = [0.5 1 0]; %[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
c_sm = [0 0.5 1]; % ;[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
suptitle('pRF center distribution')
for roi_idx = 1:num_roi
    subplot(6,2,roi_idx);
    scatter(ROI_params.prf_params{roi_idx}.x,ROI_params.prf_params{roi_idx}.y,[],c,'.');  hold on;
    scatter(ROI_params.prf_params{roi_idx}.x_sm,ROI_params.prf_params{roi_idx}.y_sm,[],c_sm,'.');
    
    tmp = strsplit(ROI_params.rois{roi_idx},'_');
    legend([tmp(2) strcat(tmp(2),'sm')],'Location','northeastoutside')
    legend('Location','NorthWestoutside')
    axis([-10 10 -10 10])
   
end

%% 3) building mrMesh and displaying parameters on the mesh

surf_visualize = 0;
if surf_visualize ==1
    % prf parameters on mrVista surface
    hvol = meshBuild(hvol,'left'); MSH = meshVisualize(viewGet(hvol,'Mesh')); hvol = viewSet(hvol, 'Mesh', MSH); clear MSH;
    % Smooth the mesh
    hvol = viewSet(hvol, 'Mesh', meshSmooth( viewGet(hvol, 'Mesh'), 1));
    
    % update map values
    prf_param = ecc;
    thr = sm_mask & ecc<20;
    
    map_val = nan(size(prf_param));
    map_val(thr) = prf_param(thr);
    
    hvol = viewSet(hvol,'displaymode','map');
    hvol = viewSet(hvol,'map',{map_val});
    hvol.ui.mapMode = setColormap(hvol.ui.mapMode,'hsvCmap');
    
    % different colormaps for phase values
    
    % Update mesh
    hvol = meshColorOverlay(hvol,1);
    
end

