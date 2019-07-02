function mprf_pRF_sm_fig(subjID, dirPth)
% Function to plot parameters on mrVista surface after mprf_pRF_sm
%                  
% INPUTS:
%   subjID      :   subject name (string)
%   dirPth      :   paths locating subject's data and files (struct, see loadPaths.m)
%  
% plots after mprf_pRF_sm
% 1) distribution of variance explained values for all voxels
% 2) prf size vs eccentricity (different rois)
%                             (original and smoothed)
% 3) mrVista surface plots (original and smoothed)


%% ------------
%   File paths
% -------------
anat_dir    = fullfile(dirPth.fmri.mrvPth, '3DAnatomy');
coords_path = fullfile(dirPth.fmri.mrvPth, 'Gray');
roi_path    = fullfile(anat_dir, 'ROIs');

% Directories to load results
prf_dir_mrv     = dirPth.fmri.saveDataPth_prfMrv;
prf_data_mrVmat = fullfile(prf_dir_mrv,'mat');      % original and smoothed pRF parameters in mrVista space in .mat
prf_data_mrVNif = fullfile(prf_dir_mrv, 'nifti');   % original and smoothed pRF parameters in mrVista space in .nii

data_file = dir(fullfile(prf_data_mrVmat,'*.mat'));
load(fullfile(prf_data_mrVmat,data_file.name));


% Other settings to define
surf_visualize = 1;
Var_Exp_Thr    = 0.4;
Ecc_Thr        = [0 10]; % degrees

%% --------------------
%   Get prf parameters
%  --------------------

% original prf parameters
ve          = prf_par_exp.varexplained;
sigma       = prf_par_exp.sigma;
x           = prf_par_exp.x;
y           = prf_par_exp.y;

% Convert to polar angle and eccentricity
[pol,ecc]   = cart2pol(x,y);
pol         = mod(pol,2*pi);
beta        = prf_par_exp.beta;

% Smoothed prf parameters
sigma_sm        = prf_par_exp.sigma_smoothed;
x_sm            = prf_par_exp.x_smoothed; y_sm = prf_par_exp.y_smoothed;
recomp_beta     = prf_par_exp.recomp_beta;

% Again, convert to polar angle and eccentricity
[pol_sm,ecc_sm] = cart2pol(x_sm,y_sm);
pol_sm          = mod(pol_sm, 2*pi);

%% ----------------------
%   Visualize some stats
%  ----------------------

% Variance explained
figure(1); set(gcf, 'Color', 'w', 'Position', [102   999   1920   400]); box off;
hist(ve,100);
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
title('Variance explained of all pRF model fits')
xlabel('Variance explained'); ylabel('Frequency')
set(gca, 'FontSize', 15, 'TickDir', 'out')

% prf size vs eccentricity (all voxels)
figure(2),set(gcf, 'Color', 'w', 'Position', [10   10   1920/3   1080/2])
sm_mask = ve > 0;
% ve_thr_mask = ve > 0.4;
c = sm_mask;
scatter(ecc,sigma,[],c);
title('voxels used for smoothing (ve>0)');
xlabel('eccentricity');
ylabel('prf size');
set(gca, 'FontSize', 15, 'TickDir', 'out')
legend('Included'); legend boxoff

%% ------------------------------------------------
%   Separate prf params for different ROIs
%  ------------------------------------------------

% Load Wang Atlas ROIs compatible with mrVista
if strcmp(subjID, 'wlsubj004')
    rois = {'WangAtlas_V1_Combined','WangAtlas_V2_Combined','WangAtlas_V3_Combined',...
        'WangAtlas_IPS_Combined','WangAtlas_LO_Combined','WangAtlas_V3AB_Combined',...
        'WangAtlas_MT_Combined','WangAtlas_VO1','WangAtlas_VO2','WangAtlas_FEF',...
        'WangAtlas_hV4'};
else
    d = dir(fullfile(roi_path, 'WangAtlas_*'));
    for rr = 1:length(d)
        thisRoi = strsplit(d(rr).name, '.');
        rois{rr} = thisRoi{1};
    end
end

rois = rois';

num_roi   = length(rois);
roi_fname = cell(num_roi,1);

for roi_idx = 1:num_roi
    roi_fname{roi_idx,1} = fullfile(roi_path,[rois{roi_idx} '.mat']);
end

ROI_params = table(rois,roi_fname);

% Load coordinate file
coordsFile = fullfile(coords_path,'coords.mat');
load(coordsFile);

% Apply these thresholds on the pRF parameters for both the conditions
model_data_thr = cell(num_roi,1);

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


% Plot pRF size vs ecc
figure(4),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
c = [0.5 1 0]; %[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
c_sm = [0 0.5 1]; % ;[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
title('pRF size vs ecc')
for roi_idx = 1:num_roi
    subplot(num_roi/2,2,roi_idx);
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
title('Sigma')
ncols = 2;
nrows = ceil(num_roi/ncols);

for roi_idx = 1:num_roi
    subplot(nrows,ncols,roi_idx);
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
title('Sigma smoothed')
for roi_idx = 1:num_roi
    subplot(nrows,ncols,roi_idx);
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
title('Beta')
for roi_idx = 1:num_roi
    subplot(nrows,ncols,roi_idx);
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
title('Recomputed beta')
for roi_idx = 1:num_roi
    subplot(nrows,ncols,roi_idx);
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
title('pRF center distribution')
for roi_idx = 1:num_roi
    subplot(nrows,ncols,roi_idx);
    scatter(ROI_params.prf_params{roi_idx}.x,ROI_params.prf_params{roi_idx}.y,[],c,'.');  hold on;
    scatter(ROI_params.prf_params{roi_idx}.x_sm,ROI_params.prf_params{roi_idx}.y_sm,[],c_sm,'.');
    
    tmp = strsplit(ROI_params.rois{roi_idx},'_');
    legend([tmp(2) strcat(tmp(2),'sm')],'Location','northeastoutside')
    legend('Location','NorthWestoutside')
    axis([-10 10 -10 10])
    
end

%% --------------------------------------------------------------------
%   Load mrMesh and display Wang Rois and prf parameters on the mesh
%  --------------------------------------------------------------------

if surf_visualize
    
    % Go to vista session and open a mrVista Gray window
    cd(dirPth.fmri.mrvPth)
    vw = mrVista('3');
    
    % Load rh and lh mesh
    mesh1 = fullfile(anat_dir, 'Left', '3DMeshes', 'Left_inflated.mat');
    mesh2 = fullfile(anat_dir, 'Right', '3DMeshes', 'Right_inflated.mat');
    [vw, OK] = meshLoad(vw, mesh1, 1); if ~OK, error('Mesh server failure'); end
    [vw, OK] = meshLoad(vw, mesh2, 1); if ~OK, error('Mesh server failure'); end
    
    % Get prf parameters saved as nifti's
    prfParams = {'eccentricity', 'eccentricity_smoothed', 'polar_angle', 'polar_angle_smoothed', ...
        'sigma', 'sigma_smoothed', 'varexplained', 'beta', 'recomp_beta'};
    
    % Load and draw Wang ROIs
    for idx = 1:num_roi
        roiFile = sprintf('%s',rois{idx});
        vw = loadROI(vw, roiFile);
        sprintf('(%s): Loaded ROI: %s \n', mfilename, roiFile)
    end
    
    vw = viewSet(vw, 'ROI draw method', 'perimeter');
    vw = refreshScreen(vw);
    vw = meshUpdateAll(vw);
    
    % Load and draw prf parameters
    for ii = 1:length(prfParams)
        
        sprintf('(%s):  Visualizing %s on mrVista surface \n', mfilename, prfParams{ii})
        
        prfParamNifti =  fullfile(prf_data_mrVNif, sprintf('%s.nii.gz', prfParams{ii}));
        
        vw = viewSet(vw,'displaymode','map');
        vw = loadParameterMap(vw,prfParamNifti);
        
        % Add smoothing/ve mask?
        vw.map{1} = vw.map{1}.*sm_mask;
        
        % Set colormap and limits
        vw.ui.mapMode = setColormap(vw.ui.mapMode, 'hsvTbCmap');
        
        switch prfParam
            case {'polar_angle', 'polar_angle_smoothed'}
                vw = viewSet(vw, 'mapwin', [eps 180]);
                vw = viewSet(vw, 'mapclip', [eps 180]);
                vw.ui.mapMode = setColormap(vw.ui.mapMode, 'hsvCmap');
            case {'eccentricity', 'eccentricity_smoothed'}
                vw = viewSet(vw, 'mapwin', [eps 20]);
                vw = viewSet(vw, 'mapclip', [eps 20]);
            case {'beta', 'recomp_beta'}
                vw = viewSet(vw, 'mapwin', [eps 20]);
                vw = viewSet(vw, 'mapclip', [eps 20]);
        end
        
        % Update views
        vw = refreshScreen(vw);
        vw = meshUpdateAll(vw);       
        
        % Copy the mesh to a Matlab figure
        % fH = figure('Color', 'w');
        % imagesc(mrmGet(viewGet(vw, 'Mesh'), 'screenshot')/255); axis image; axis off;
    
    end
end


end

