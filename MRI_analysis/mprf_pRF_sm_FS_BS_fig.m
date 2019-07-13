function mprf_pRF_sm_FS_BS_fig(dirPth,opt)
% plots after mprf_pRF_sm
% 1) distribution of variance explained values for all voxels
% 2) prf size vs eccentricity (different rois)
%                             (original and smoothed)
% 3) surface plots (original and smoothed)

% File paths
% ----------
prf_dir_BS = dirPth.fmri.saveDataPth_prfBS;
roi_dir_BS = dirPth.fmri.saveDataPth_roiBS;
% ----------

% Get directory to save images
saveDir = fullfile(dirPth.fmri.saveDataPth_prfBS, 'figs');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

% variance explained threshold for selecting voxels 
Var_Exp_Thr = opt.varExplThresh(1);

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
fH1 = figure(201); set(gcf, 'Color', 'w', 'Position', [102   999   1920   400])
hist(bs_prf.varexplained{1},100);
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
title('variance explained of prf fits');
xlabel('variance explained');
print(fH1, fullfile(saveDir,'variance_explained'), '-dpng');

% prf size vs eccentricity
%-------------------------
% All voxels
fH2 = figure(202); set(gcf, 'Color', 'w', 'Position', [10   10   1920/3   1080/2])
sm_mask = bs_prf.varexplained{1} > 0;
scatter(bs_prf.eccentricity{1},bs_prf.sigma{1},[],[0.5 0.5 0.5]); hold on;
scatter(bs_prf.eccentricity{1}(sm_mask),bs_prf.sigma{1}(sm_mask),[],[0.5 1 0.5]); hold on;
title('voxels used for smoothing (ve>0)');
xlabel('eccentricity');
ylabel('prf size');
print(fH2, fullfile(saveDir,'voxels used for smoothing'), '-dpng');

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
idx_thr = bs_prf.varexplained{1} > Var_Exp_Thr;% & bs_prf.eccentricity{1} < Ecc_Thr(2);


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
fH3 = figure(203); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
c = [0.5 1 0]; %[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
c_sm = [0 0.5 1]; % ;[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
title('pRF size vs ecc')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
  
    scatter(bs_prf_roi.eccentricity{roi_idx},bs_prf_roi.sigma{roi_idx},[],c,'*');  hold on;
    scatter(bs_prf_roi.eccentricity_smoothed{roi_idx},bs_prf_roi.sigma_smoothed{roi_idx},[],c_sm,'*');
    
    tmp = bs_prf_roi.ROIs{roi_idx};
    legend([{tmp} {strcat(tmp,'sm')}],'Location','NorthWestOutside')
    ylim([0 15]);
end
print(fH3, fullfile(saveDir,'pRF_size_eccentricity'), '-dpng');

% Sigma histogram
%-----------------------
fH4 = figure(204); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
title('Sigma')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(bs_prf_roi.sigma{roi_idx},100)
    tmp = bs_prf_roi.ROIs{roi_idx};
    title(tmp)   
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end
print(fH4, fullfile(saveDir,'sigma'), '-dpng');

% Smoothed sigma histogram
%-----------------------
fH5 = figure(205); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
title('Sigma smoothed')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(bs_prf_roi.sigma_smoothed{roi_idx},100)
    tmp = bs_prf_roi.ROIs{roi_idx};
    title(tmp)   
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end
print(fH5, fullfile(saveDir,'sigma_smoothed'), '-dpng');

% Beta histogram
%-----------------------
fH6 = figure(206); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
title('Beta')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(bs_prf_roi.beta{roi_idx},100)
    tmp = bs_prf_roi.ROIs{roi_idx};
    title(tmp)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end
print(fH6, fullfile(saveDir,'beta'), '-dpng');

% recomputed beta histogram
%-----------------------
fH7 = figure(207); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
title('Recomputed beta')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    hist(bs_prf_roi.recomp_beta{roi_idx},100)
    tmp = bs_prf_roi.ROIs{roi_idx};
    title(tmp)  
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end
print(fH7, fullfile(saveDir,'recomp_beta'), '-dpng');

fH8 = figure(208); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
c = [0.5 1 0];
c_sm = [0 0.5 1];
title('pRF center distribution')
for roi_idx = 1:num_roi
    subplot(num_row,num_col,roi_idx);
    scatter(bs_prf_roi.x{roi_idx},bs_prf_roi.y{roi_idx},[],c,'.');  hold on;
    scatter(bs_prf_roi.x_smoothed{roi_idx},bs_prf_roi.y_smoothed{roi_idx},[],c_sm,'.');
    %axis image
    xlim([-10 10])
    ylim([-10 10])
    
    tmp = bs_prf_roi.ROIs{roi_idx};
    legend([{tmp} {strcat(tmp,'sm')}],'Location','northeastoutside')
    legend('Location','NorthWestoutside');
end
print(fH8, fullfile(saveDir,'prf_center_distribution'), '-dpng');

%% 3) building mrMesh and displaying parameters on the mesh

if opt.surfVisualize ==1
    close all;
   
    % Get directory to save images
    saveDir = fullfile(dirPth.fmri.saveDataPth_prfBS, 'figs');
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
   
    %-----------------------------------------------
    % Visualize pRF parameters on brainstorm surface
    %-----------------------------------------------
    % Load rh and lh freesurfer surface files
    surfaces_to_load = {'pial'};
          
    for idx_surf = 1:length(surfaces_to_load)
        % left hemisphere
        cur_surf = surfaces_to_load{idx_surf};
        mprf_VisualizeDataOnBrainstormSurface(dirPth,cur_surf,saveDir);
        
        vis_roi=0;
        if vis_roi==1
            % Hide rois in gray view when loading
            vw = viewSet(vw, 'hide gray rois', true);
            
            % Load and draw Wang ROIs
            for idx = 1:num_roi
                roiFile = sprintf('%s',rois{idx});
                vw = loadROI(vw, roiFile);
                fprintf('(%s): Loaded ROI: %s \n', mfilename, roiFile)
            end
            
            vw = viewSet(vw, 'ROI draw method', 'perimeter');
            vw = refreshScreen(vw);
            vw = meshUpdateAll(vw);
        end
        
    end
    
    
end


