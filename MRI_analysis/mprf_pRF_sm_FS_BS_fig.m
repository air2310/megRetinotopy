% plots after mprf_pRF_sm
% 1) distribution of variance explained values for all voxels
% 2) prf size vs eccentricity (different rois)
%                             (original and smoothed)
% 3) surface plots (original and smoothed)

%% 1) Variance explained distribution
subjid = 'wlsubj004';

% directories required
Anat_dir = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Data/Anatomy/%s',subjid);
bs_prf_data = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Quality_check/%s/prf_data/surface/brainstorm',subjid);
bs_roi_dir = sprintf('/mnt/storage_2/projects/MEG/Retinotopy/Quality_check/%s/rois/surface/brainstorm',subjid);

% read prf parameters (original and smoothed) in brainstorm space
% - variance explained - 
% - sigma
% - x & y
% - beta & recomp_beta
% - polar angle

% ROIs

surfaces_to_load = {'pial'};
pname = bs_prf_data;
w_export = 'roi';

prf_files = dir(fullfile(pname,strcat(surfaces_to_load{1},'.*')));
bs_prf = table(surfaces_to_load);
% Loop over the LHs files (i.e. all the parameters on the lh surfaces):
for nn = 1:length(prf_files)
    
    cur_prf_file = prf_files(nn).name;
    par_name = strsplit(cur_prf_file,'.');
    
    % load and concatenate:
    prf_data = read_curv(fullfile(pname,cur_prf_file));
    
    add_t_1 = table({prf_data},'variableNames',par_name(2));clos
    
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
pname = bs_roi_dir;

w_export = 'roi';


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



