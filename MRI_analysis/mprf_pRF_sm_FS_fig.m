function mprf_pRF_sm_FS_fig(dirPth)

% plots after mprf_pRF_sm
% 1) distribution of variance explained values for all voxels
% 2) prf size vs eccentricity (different rois)
%                             (original and smoothed)
% 3) surface plots (original and smoothed)

% File paths
% ----------
rootDir = dirPth.rootPth;
prf_dir_FS = strcat(rootDir,dirPth.fmri.saveDataPth_prfFS(2:end));
fs_roi_dir = strcat(rootDir,dirPth.fmri.saveDataPth_roiFS(2:end));
% ----------

pname = prf_dir_FS;

surfaces_to_load = {'pial'};
hs_to_load = {'lh','rh'};

lh_files = dir(fullfile(pname,strcat(hs_to_load{1},'.*')));
both_data = cell(length(lh_files),1);
fs_prf = table(surfaces_to_load);
for nn = 1:length(lh_files)
    cur_lh_file = lh_files(nn).name;
    [~,~,par_name] = fileparts(cur_lh_file);
    
    par_name = par_name(2:end);
    % Find the corresponding rh file:
    cur_rh_file = ['rh.' par_name];
    
    % load and concatenate:
    both_data = [read_curv(fullfile(pname,cur_lh_file));...
        read_curv(fullfile(pname,cur_rh_file))];
    
    add_t_1 = table({both_data},'variableNames',{par_name});
    
    fs_prf = [fs_prf add_t_1];   
end


% histogram of variance explained 
%--------------------------------
figure(101); set(gcf, 'Color', 'w', 'Position', [102   999   1920   400])
hist(fs_prf.varexplained{1},100);
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
title('variance explained of prf fits');
xlabel('variance explained');

% prf size vs eccentricity
%-------------------------
% All voxels
figure(102),set(gcf, 'Color', 'w', 'Position', [10   10   1920/3   1080/2])
sm_mask = fs_prf.varexplained{1} > 0;
c = sm_mask;
scatter(fs_prf.eccentricity{1},fs_prf.sigma{1},[],c);
title('voxels used for smoothing (ve>0)');
xlabel('eccentricity');
ylabel('prf size');


% Check how ROIs are loaded
pname = fs_roi_dir;

% ROIs
lh_files = dir(fullfile(pname,'lh.*'));
fs_roi = table(surfaces_to_load);
% Loop over the LHs files (i.e. all the parameters on the lh surfaces):
for nn = 1:length(lh_files)
    
    cur_lh_file = lh_files(nn).name;
    [~,~,par_name] = fileparts(cur_lh_file);
    par_name = par_name(2:end);
    % Find the corresponding rh file:
    cur_rh_file = ['rh.' par_name];
    
    % load and concatenate:
    both_data = [read_curv(fullfile(pname,cur_lh_file));...
        read_curv(fullfile(pname,cur_rh_file))];
    
    add_t_1 = table({both_data},'variableNames',{par_name});
    
    fs_roi = [fs_roi add_t_1];
end
        

% fits for different visual areas
Var_Exp_Thr = 0.4;
Ecc_Thr = [0 10];
idx_thr = fs_prf.varexplained{1} > Var_Exp_Thr & fs_prf.eccentricity{1} < Ecc_Thr(2);


fs_prf_roi = table(fs_roi.Properties.VariableNames(2:end)','VariableNames',{'ROIs'});
for idx_param = 1:(size(fs_prf,2)-1)
    cur_param = nan(size(fs_prf{1,idx_param+1}{1}));
    cur_param(idx_thr) = fs_prf{1,idx_param+1}{1}(idx_thr);
    cur_param_name = fs_prf(1,idx_param+1); 
    clear add_t_2       
    for idx_roi = 1:(size(fs_roi,2)-1)
 
        cur_roi = find(~isnan(fs_roi{1,idx_roi+1}{1}));    
        cur_roi_name = fs_roi(1,idx_roi+1); 
        
        cur_param_val_thr(idx_roi,1) = {cur_param(cur_roi)}; 
        
        
    end
    add_t_2 = table(cur_param_val_thr,'VariableNames',cur_param_name.Properties.VariableNames);
    
    fs_prf_roi = [fs_prf_roi add_t_2];
end

num_roi = length(fs_prf_roi.ROIs);
num_row = 5;
num_col = 3;

% pRF size vs eccentricity
%----------------------
figure(103),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
c = [0.5 1 0]; %[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
c_sm = [0 0.5 1]; % ;[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
suptitle('pRF size vs ecc')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
  
    scatter(fs_prf_roi.eccentricity{roi_idx},fs_prf_roi.sigma{roi_idx},[],c,'*');  hold on;
    scatter(fs_prf_roi.eccentricity_smoothed{roi_idx},fs_prf_roi.sigma_smoothed{roi_idx},[],c_sm,'*');
    
    tmp = fs_prf_roi.ROIs{roi_idx};
    legend([{tmp} {strcat(tmp,'sm')}],'Location','NorthWestOutside')
    ylim([0 15]);
end

% Sigma histogram
%-----------------------
figure(104),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Sigma')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(fs_prf_roi.sigma{roi_idx},100)
    tmp = fs_prf_roi.ROIs{roi_idx};
    title(tmp)   
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end

% Smoothed sigma histogram
%-----------------------
figure(105),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Sigma smoothed')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(fs_prf_roi.sigma_smoothed{roi_idx},100)
    tmp = fs_prf_roi.ROIs{roi_idx};
    title(tmp)   
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end


% Beta histogram
%-----------------------
figure(106),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Beta')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(fs_prf_roi.beta{roi_idx},100)
    tmp = fs_prf_roi.ROIs{roi_idx};
    title(tmp)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end

% recomputed beta histogram
%-----------------------
figure(107),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
suptitle('Recomputed beta')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(fs_prf_roi.recomp_beta{roi_idx},100)
    tmp = fs_prf_roi.ROIs{roi_idx};
    title(tmp)  
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end

figure(108),set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
c = [0.5 1 0]; %[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
c_sm = [0 0.5 1]; % ;[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
suptitle('pRF center distribution')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    scatter(fs_prf_roi.x{roi_idx},fs_prf_roi.y{roi_idx},[],c,'.');  hold on;
    scatter(fs_prf_roi.x_smoothed{roi_idx},fs_prf_roi.y_smoothed{roi_idx},[],c_sm,'.');
    %axis image
    xlim([-10 10])
    ylim([-10 10])
    
    tmp = fs_prf_roi.ROIs{roi_idx};
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

