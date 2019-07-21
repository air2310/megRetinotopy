function mprf_pRF_sm_FS_fig(dirPth,opt)

% plots after mprf_pRF_sm
% 1) distribution of variance explained values for all voxels
% 2) prf size vs eccentricity (different rois)
%                             (original and smoothed)
% 3) surface plots (original and smoothed)

%% ------------
%   File paths
% -------------

prfFSDir = dirPth.fmri.saveDataPth_prfFS;
roiFSDir = dirPth.fmri.saveDataPth_roiFS;
% ----------

% Get directory to save images
saveDir = fullfile(dirPth.fmri.saveDataPth_prfFS, 'figs');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

% load left hemi files
lhFiles = dir(fullfile(prfFSDir,'lh.*'));

prfDataFS = struct();

for nn = 1:length(lhFiles)
    
    curLHFile = lhFiles(nn).name;
    tmp = str_split(curLHFile, '.');
    prfName = tmp{2};
    
    % Find the corresponding rh file:
    curRHFile = ['rh.' prfName];
    
    % load and concatenate:
    bothFSHemiData = [read_curv(fullfile(prfFSDir,curLHFile));...
        read_curv(fullfile(prfFSDir,curRHFile))];
    
    prfDataFS.(prfName) = bothFSHemiData;
    
end


% histogram of variance explained
%--------------------------------
fH1 = figure(101); set(gcf, 'Color', 'w', 'Position', [102   999   1920   400])
hist(prfDataFS.varexplained,100);
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
title('variance explained of prf fits on FS surface');
xlabel('variance explained');
print(fH1, fullfile(saveDir,'variance_explained'), '-dpng');

% prf size vs eccentricity
%-------------------------
% All voxels
fH2 = figure(102);set(gcf, 'Color', 'w', 'Position', [10   10   1920/3   1080/2])
sm_mask = prfDataFS.varexplained > 0;
c = sm_mask;
scatter(prfDataFS.eccentricity,prfDataFS.sigma,[],c);
title('All voxels with ve >0');
xlabel('eccentricity');
ylabel('prf size');
print(fH2, fullfile(saveDir,'voxels_vethresh'), '-dpng');


%% Load FS ROIs

lhROIFiles = dir(fullfile(roiFSDir,'lh.*'));
roiDataFS = struct();

% Loop over the LH roi files
for nn = 1:length(lhROIFiles)
    
    curLHFile = lhROIFiles(nn).name;
    tmp = str_split(curLHFile, '.');
    roiName = tmp{2};
    
    % Find the corresponding rh file:
    curRHFile = ['rh.' roiName];
    
    % load and concatenate:
    bothHemiROI = [read_curv(fullfile(roiFSDir,curLHFile)); ...
        read_curv(fullfile(roiFSDir,curRHFile))];
    
    if strcmp(roiName,'wang2015_atlas')
        wangROINames = {'V1v', 'V1d', 'V2v', 'V2d', 'V3v', 'V3d', 'hV4', 'VO1', 'VO2', 'PHC1', 'PHC2', ...
            'TO2', 'TO1', 'LO2', 'LO1', 'V3B', 'V3A', 'IPS0', 'IPS1', 'IPS2', 'IPS3', 'IPS4', ...
            'IPS5', 'SPL1', 'FEF'};
        
        roisToCombine = {'V1','V2','V3','IPS','LO','TO','VO', 'PHC'};
        roisIndiv     = find(~contains(wangROINames,roisToCombine));
        
        % Combine data of ROIs
        for rc = 1:length(roisToCombine)
            match = find(contains(wangROINames,roisToCombine(rc)));
            
            tmpData = [];
            for ii = 1:length(match)
                indices = find(bothHemiROI==match(ii));
                tmpData = [tmpData; indices];
            end
            
            roiDataFS.(roisToCombine{rc}) = unique(tmpData, 'rows');
        end
        % Add ROIs that do not need to be combined
        for rc = roisIndiv
            roiDataFS.(wangROINames{rc}) =  bothHemiROI(bothHemiROI==rc);
        end 
    else
        roiDataFS.(roiName) = find(bothHemiROI);
    end
end



% Sample prf data for different visual areas
vemask     = (prfDataFS.varexplained > opt.varExplThresh(1) & prfDataFS.varexplained < opt.varExplThresh(2));
eccenmask  = (prfDataFS.eccentricity > opt.eccThresh(1) & prfDataFS.eccentricity < opt.eccThresh(2));
prfROIData = struct();
fnData     = fieldnames(prfDataFS);
fnRoi      = fieldnames(roiDataFS);

for p = 1:length(fnData)
    
    allroitmp = [];
    
    for r = 1:length(fnRoi)
        
        thisFieldName = fnData{p};
        data = prfDataFS.(fnData{p});
        
        if strcmp(thisFieldName,'varexplained')
            data = data(roiDataFS.(fnRoi{r}));
            prfROIData.(fnRoi{r}).(fnData{p}) = data(data > opt.varExplThresh(1) & data < opt.varExplThresh(2));
            
        elseif (strcmp(thisFieldName,'beta') || strcmp(thisFieldName,'recomp_beta'))
            data = data(roiDataFS.(fnRoi{r}));
            thresh = prctile(data, opt.betaPrctileThresh);
            betamask = ((data > thresh(1)) & (data < thresh(2)));
            prfROIData.(fnRoi{r}).(fnData{p}) = data(betamask);
        else
            
            % Mask
            data = data(roiDataFS.(fnRoi{r}));
            data = data((vemask(roiDataFS.(fnRoi{r})) & eccenmask(roiDataFS.(fnRoi{r}))));
            
            prfROIData.(fnRoi{r}).(fnData{p}) = data;
        end
        allroitmp  = [allroitmp; data];
    end
    prfROIData.allROI.(thisFieldName) = allroitmp;
end

fnRoi(end+1) = {'allROI'};
numRois = length(fnRoi);


nRow = 5;
nCol = 3;

% pRF size vs eccentricity
%----------------------
fH3 = figure(103);set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
c = [0.5 1 0]; %[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
c_sm = [0 0.5 1]; % ;[0.5 0 0];[1 0 0];[0 0.5 0];[0 1 0];[0 0 0.5];[0 0 1];[0.5 0.5 0.5];[1 0.5 0.5];[0.5 1 0.5];[0.5 0.5 1]];
title('pRF size vs ecc')
for roi_idx = 1:numRois
    subplot(nRow,nCol,roi_idx);
    
    scatter(prfROIData.(fnRoi{roi_idx}).eccentricity, prfROIData.(fnRoi{roi_idx}).sigma,[],c,'*');  hold on;
    scatter(prfROIData.(fnRoi{roi_idx}).eccentricity_smoothed,prfROIData.(fnRoi{roi_idx}).sigma_smoothed,[],c_sm,'*');
    
    title(fnRoi{roi_idx}, 'FontSize', 15)
    legend('original','smoothed','Location','NorthWest')
    legend('Location','NorthWest')
    ylim([0 opt.eccThresh(2)]); xlim([0 opt.eccThresh(2)]); axis square
    ylabel('PRF sigma (deg)'); xlabel('PRF eccen (deg');
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH3, fullfile(saveDir,'pRF_size_eccentricity'), '-dpng');

% Sigma histogram
%-----------------------
fH4 = figure(104); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);

for roi_idx = 1:numRois
    subplot(nRow,nCol,roi_idx);
    hist(prfROIData.(fnRoi{roi_idx}).sigma,100)
    
    title(fnRoi{roi_idx}, 'FontSize', 15)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    xlabel('Sigma'); ylabel('Frequency')
end
print(fH4, fullfile(saveDir,'sigma'), '-dpng');

% Smoothed sigma histogram
%-----------------------
fH5 = figure(105); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);

for roi_idx = 1:numRois
    subplot(nRow,nCol,roi_idx);
    hist(prfROIData.(fnRoi{roi_idx}).sigma_smoothed,100)
    
    title(fnRoi{roi_idx}, 'FontSize', 15)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    xlabel('Sigma smoothed'); ylabel('Frequency')
end
print(fH5, fullfile(saveDir,'sigma_smoothed'), '-dpng');

% Beta histogram
%-----------------------
fH6 = figure(106); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);

for roi_idx = 1:numRois
    subplot(nRow,nCol,roi_idx);
    hist(prfROIData.(fnRoi{roi_idx}).beta,100)
    
    title(fnRoi{roi_idx}, 'FontSize', 15)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    xlabel('Beta'); ylabel('Frequency')
end
print(fH6, fullfile(saveDir,'beta'), '-dpng');

% recomputed beta histogram
%-----------------------
fH7 = figure(107); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
for roi_idx = 1:numRois
    subplot(nRow,nCol,roi_idx);
    hist(prfROIData.(fnRoi{roi_idx}).recomp_beta,100)
    title(fnRoi{roi_idx}, 'FontSize', 15)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    xlabel('Recomp Beta'); ylabel('Frequency')
end
print(fH7, fullfile(saveDir,'recomp_beta'), '-dpng');

fH8 = figure(108); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
c = [0.5 1 0]; 
c_sm = [0 0.5 1]; 
title('pRF center distribution')
for roi_idx = 1:numRois
    subplot(nRow,nCol,roi_idx);
    scatter(prfROIData.(fnRoi{roi_idx}).x,prfROIData.(fnRoi{roi_idx}).y,[],c,'.');  hold on;
    scatter(prfROIData.(fnRoi{roi_idx}).x_smoothed,prfROIData.(fnRoi{roi_idx}).y_smoothed,[],c_sm,'.');
    %axis image
    xlim([-10 10])
    ylim([-10 10])
    
    title(fnRoi{roi_idx}, 'FontSize', 15)
    legend({'Orig', 'Smoothed'}, 'Location','NorthWestoutside')
    
    
end
print(fH8, fullfile(saveDir,'prf_center_distribution'), '-dpng');

%% --------------------------------------------------------------------
%   Load mrMesh and display Wang Rois and prf parameters on the mesh
%  --------------------------------------------------------------------

if opt.surfVisualize
    close all;
    
    % Get directory to save images
    saveDir = fullfile(dirPth.fmri.saveDataPth_prfFS, 'figs');
    saveDir = fullfile(dirPth.fmri.saveDataPth_prfFS, 'figs','rh');
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    
    %-----------------------------------------------
    % Visualize pRF parameters on freesurfer surface
    %-----------------------------------------------
    % Load rh and lh freesurfer surface files
    hs_to_load = {'lh','rh'};
    surfaces_to_load = {'pial'};
    
     for cur_surf = 1:length(surfaces_to_load)
        for cur_hs = 1:length(hs_to_load)
            
            cur_surf_to_load = strcat(hs_to_load{cur_hs},'.',surfaces_to_load{cur_surf});
            cur_roi_to_load = strcat(hs_to_load{cur_hs},'.',surfaces_to_load{cur_surf});
            
            mprf_VisualizeDataOnFreesurferSurface(dirPth, cur_surf_to_load, saveDir);
            
            mprf_VisualizeRoiOnFreesurferSurface(dirPth, cur_surf_to_load,saveDir);
        end
    end
end

end