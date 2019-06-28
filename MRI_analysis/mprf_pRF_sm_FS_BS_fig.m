function mprf_pRF_sm_FS_BS_fig(dirPth)
% plots after mprf_pRF_sm
% 1) distribution of variance explained values for all voxels
% 2) prf size vs eccentricity (different rois)
%                             (original and smoothed)
% 3) surface plots (original and smoothed)

% File paths
% ----------
rootDir = dirPth.rootPth;
prf_dir_BS = strcat(rootDir,dirPth.fmri.saveDataPth_prfBS(2:end));
roi_dir_BS = strcat(rootDir,dirPth.fmri.saveDataPth_roiBS(2:end));
% ----------

% ROIs
surfaces_to_load = {'pial'};
pname = prf_dir_BS;

prf_files = dir(fullfile(pname,strcat(surfaces_to_load{1},'.*')));
bs_prf = table(surfaces_to_load);
% Loop over the LHs files (i.e. all the parameters on the lh surfaces):
for nn = 1:length(prf_files)
    
    cur_prf_file = prf_files(nn).name;
    par_name = strsplit(cur_prf_file,'.');
    
    % load and concatenate:
    prf_data = read_curv(fullfile(pname,cur_prf_file));
    
    add_t_1 = table({prf_data},'variableNames',par_name(2));
    
    bs_prf = [bs_prf add_t_1];
end


% histogram of variance explained 
%--------------------------------
figure(201); set(gcf, 'Color', 'w', 'Position', [102   999   1920   400])
hist(bs_prf.varexplained{1},100);
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
title('variance explained of prf fits');
xlabel('variance explained');

% prf size vs eccentricity
%-------------------------
% All voxels
figure(202),set(gcf, 'Color', 'w', 'Position', [10   10   1920/3   1080/2])
sm_mask = bs_prf.varexplained{1} > 0;
scatter(bs_prf.eccentricity{1},bs_prf.sigma{1},[],[0.5 0.5 0.5]); hold on;
scatter(bs_prf.eccentricity{1}(sm_mask),bs_prf.sigma{1}(sm_mask),[],[0.5 1 0.5]); hold on;
title('voxels used for smoothing (ve>0)');
xlabel('eccentricity');
ylabel('prf size');


% Check how ROIs are loaded
pname = roi_dir_BS;

% ROIs
roi_files = dir(fullfile(pname,strcat(surfaces_to_load{1},'.*')));
bs_roi = table(surfaces_to_load);
% Loop over the LHs files (i.e. all the parameters on the lh surfaces):
for nn = 1:length(roi_files)
    
    cur_roi_file = roi_files(nn).name;
    par_name = strsplit(cur_roi_file,'.');
    
    % load and concatenate:
    roi_data = read_curv(fullfile(pname,cur_roi_file));
    
    add_t_1 = table({roi_data},'variableNames',par_name(2));
    
    bs_roi = [bs_roi add_t_1];
end
        

% fits for different visual areas
Var_Exp_Thr = 0.4;
Ecc_Thr = [0 10];
idx_thr = bs_prf.varexplained{1} > Var_Exp_Thr & bs_prf.eccentricity{1} < Ecc_Thr(2);


bs_prf_roi = table(bs_roi.Properties.VariableNames(2:end)','VariableNames',{'ROIs'});
for idx_param = 1:(size(bs_prf,2)-1)
    cur_param = nan(size(bs_prf{1,idx_param+1}{1}));
    cur_param(idx_thr) = bs_prf{1,idx_param+1}{1}(idx_thr);
    cur_param_name = bs_prf(1,idx_param+1); 
    clear add_t_2       
    for idx_roi = 1:(size(bs_roi,2)-1)
 
        cur_roi = find(~isnan(bs_roi{1,idx_roi+1}{1}));    
        cur_roi_name = bs_roi(1,idx_roi+1); 
        
        cur_param_val_thr(idx_roi,1) = {cur_param(cur_roi)}; 
        
        
    end
    add_t_2 = table(cur_param_val_thr,'VariableNames',cur_param_name.Properties.VariableNames);
    
    bs_prf_roi = [bs_prf_roi add_t_2];
end

num_roi = length(bs_prf_roi.ROIs);
num_row = 5;
num_col = 3;

% pRF size vs eccentricity
%----------------------
figure(203),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
c = [0.5 1 0]; %[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
c_sm = [0 0.5 1]; % ;[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
suptitle('pRF size vs ecc')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
  
    scatter(bs_prf_roi.eccentricity{roi_idx},bs_prf_roi.sigma{roi_idx},[],c,'*');  hold on;
    scatter(bs_prf_roi.eccentricity_smoothed{roi_idx},bs_prf_roi.sigma_smoothed{roi_idx},[],c_sm,'*');
    
    tmp = bs_prf_roi.ROIs{roi_idx};
    legend([{tmp} {strcat(tmp,'sm')}],'Location','NorthWestOutside')
    ylim([0 15]);
end

% Sigma histogram
%-----------------------
figure(204),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Sigma')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(bs_prf_roi.sigma{roi_idx},100)
    tmp = bs_prf_roi.ROIs{roi_idx};
    title(tmp)   
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end

% Smoothed sigma histogram
%-----------------------
figure(205),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Sigma smoothed')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(bs_prf_roi.sigma_smoothed{roi_idx},100)
    tmp = bs_prf_roi.ROIs{roi_idx};
    title(tmp)   
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end


% Beta histogram
%-----------------------
figure(206),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Beta')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(bs_prf_roi.beta{roi_idx},100)
    tmp = bs_prf_roi.ROIs{roi_idx};
    title(tmp)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end

% recomputed beta histogram
%-----------------------
figure(207),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Recomputed beta')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(bs_prf_roi.recomp_beta{roi_idx},100)
    tmp = bs_prf_roi.ROIs{roi_idx};
    title(tmp)  
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end

figure(208),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
c = [0.5 1 0]; %[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
c_sm = [0 0.5 1]; % ;[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
suptitle('pRF center distribution')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    scatter(bs_prf_roi.x{roi_idx},bs_prf_roi.y{roi_idx},[],c,'.');  hold on;
    scatter(bs_prf_roi.x_smoothed{roi_idx},bs_prf_roi.y_smoothed{roi_idx},[],c_sm,'.');
    %axis image
    xlim([-10 10])
    ylim([-10 10])
    
    tmp = bs_prf_roi.ROIs{roi_idx};
    legend([{tmp} {strcat(tmp,'sm')}],'Location','northeastoutside')
    legend('Location','NorthWestoutside')
    

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


