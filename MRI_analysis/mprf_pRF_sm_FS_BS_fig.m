function mprf_pRF_sm_FS_BS_fig(dirPth,opt)
% plots after mprf_pRF_sm_FS_BS
% 1) distribution of variance explained values for all Brainstorm vertices
% 2) prf size vs eccentricity (different rois)
%                             (original and smoothed)
% 3) surface plots (original and smoothed)

% File paths
% ----------
prfBSDir = dirPth.fmri.saveDataPth_prfBS;
roiBSDir = dirPth.fmri.saveDataPth_roiBS;
% ----------

% Get directory to save images
saveDir = fullfile(dirPth.fmri.saveDataPth_prfBS, 'figs');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

% Define file names
surfacesToLoad = 'pial';
prfFiles = dir(fullfile(prfBSDir,strcat(surfacesToLoad,'.*')));

% Loop over prf pial surface files
for nn = 1:length(prfFiles)
    cur_prf_file = prfFiles(nn).name;
    paramName = strsplit(cur_prf_file,'.');
    
    if sum(cellfun(@isempty, (regexp(paramName{2}, {'mask', 'V123mask', 'wang2015_atlas'}, 'match'))))==3
        % load and concatenate:
        prfData.(paramName{2}) = read_curv(fullfile(prfBSDir,cur_prf_file));
    end
end

% Number of bins for histograms
nbins = 50;

% histogram of variance explained
%--------------------------------
fH1 = figure(201); set(gcf, 'Color', 'w', 'Position', [102   999   1920   400], 'Name', 'Variance explained of all BS vertices')
hist(prfData.varexplained,nbins);
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
title('variance explained of prf fits');
xlabel('variance explained');
print(fH1, fullfile(saveDir,'variance_explained'), '-dpng');

% prf size vs eccentricity
%-------------------------
% All voxels
fH2 = figure(202); set(gcf, 'Color', 'w', 'Position', [10, 10, 1920/3, 1080/2], 'Name', 'BS vertices with ve>0')
sm_mask = prfData.varexplained > 0;
scatter(prfData.eccentricity, prfData.sigma,[],[0.5 0.5 0.5]); hold on;
scatter(prfData.eccentricity_smoothed(sm_mask), prfData.sigma_smoothed(sm_mask),[],[0.5 1 0.5]); hold on;
title('voxels used for smoothing (ve>0)');
xlabel('Eccentricity (deg)');
ylabel('pRF size (deg)');
print(fH2, fullfile(saveDir,'voxels_used_for_smoothing'), '-dpng');


%% ROIs

% load wang rois
combineHemi = false;
combineRois = true;
roiData = loadWangROIs(roiBSDir, 'pial.wang2015_atlas', combineHemi, combineRois);

% construct 2 masks: one for Variance Explained above 10%, one with
% eccentricity below 10 deg.
vemask    = (prfData.varexplained > opt.varExplThresh(1) & prfData.varexplained < opt.varExplThresh(2));
eccenmask = (prfData.eccentricity > opt.eccThresh(1) & prfData.eccentricity < opt.eccThresh(2));

% Apply roi masks to pRF parameters
fnPrf = fieldnames(prfData);
fnRoi = fieldnames(roiData);
numRoi = length(fnRoi);

prfROIData = struct();

for p = 1:length(fnPrf)
    
    for r = 1:numRoi
        
        thisFieldName = fnPrf{p};
        data = prfData.(fnPrf{p});
        roiMask = roiData.(fnRoi{r})>0;
        
        if strcmp(thisFieldName,'varexplained')
            data = data(vemask(roiMask));
            prfROIData.(fnRoi{r}).(fnPrf{p}) = data;
            
        elseif (strcmp(thisFieldName,'beta') || strcmp(thisFieldName,'recomp_beta'))
            data = data(roiMask);
            thresh = prctile(data, opt.betaPrctileThresh);
            betamask = ((data > thresh(1)) & (data < thresh(2)));
            prfROIData.(fnRoi{r}).(fnPrf{p}) = data(betamask);
        else
            
            % Mask
            %             data = data(roiData.fnRoi{r});
            data = data((vemask(roiMask) & eccenmask(roiMask)));
            
            prfROIData.(fnRoi{r}).(fnPrf{p}) = data;
        end
    end
end



numRow = round(numRoi/5);
numCol = ceil(numRoi/numRow);

% pRF size vs eccentricity
%----------------------
fH3 = figure(203); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080],  'Name', 'BS vertices: eccen vs sigma');
c = [0.5 1 0];
c_sm = [0 0.5 1];
for roi_idx = 1:numRoi
    subplot(numRow,numCol,roi_idx);
    
    scatter(prfROIData.(fnRoi{roi_idx}).eccentricity,prfROIData.(fnRoi{roi_idx}).sigma,[],c,'*');  hold on;
    scatter(prfROIData.(fnRoi{roi_idx}).eccentricity_smoothed,prfROIData.(fnRoi{roi_idx}).sigma_smoothed,[],c_sm,'*');
    
    title(fnRoi{roi_idx});
    legend({'orig','smoothed'},'Location','NorthWestOutside')
    ylim([0 15]);
end
print(fH3, fullfile(saveDir,'pRF_size_eccentricity'), '-dpng');

% Sigma histogram
%-----------------------
fH4 = figure(204); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'BS vertices: sigma');
for roi_idx = 1:numRoi
    subplot(numRow,numCol,roi_idx);
    hist(prfROIData.(fnRoi{roi_idx}).sigma,nbins)
    title(fnRoi{roi_idx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end
print(fH4, fullfile(saveDir,'sigma'), '-dpng');

% Smoothed sigma histogram
%-----------------------
fH5 = figure(205); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'BS vertices: sigma smoothed');
for roi_idx = 1:numRoi
    subplot(numRow,numCol,roi_idx);
    hist(prfROIData.(fnRoi{roi_idx}).sigma_smoothed,nbins)
    title(fnRoi{roi_idx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end
print(fH5, fullfile(saveDir,'sigma_smoothed'), '-dpng');

% Beta histogram
%-----------------------
fH6 = figure(206); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'BS vertices:  Beta');
for roi_idx = 1:numRoi
    subplot(numRow,numCol,roi_idx);
    hist(prfROIData.(fnRoi{roi_idx}).beta,nbins)
    title(fnRoi{roi_idx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end
print(fH6, fullfile(saveDir,'beta'), '-dpng');

% recomputed beta histogram
%-----------------------
fH7 = figure(207); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'BS vertices:  Beta Recomputed');
for roi_idx = 1:numRoi
    subplot(numRow,numCol,roi_idx);
    hist(prfROIData.(fnRoi{roi_idx}).recomp_beta,nbins)
    title(fnRoi{roi_idx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
end
print(fH7, fullfile(saveDir,'recomp_beta'), '-dpng');


fH8 = figure(208); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'BS vertices:  center (x,y)');
c = [0.5 1 0];
c_sm = [0 0.5 1];
for roi_idx = 1:numRoi
    subplot(numRow,numCol,roi_idx);
    
    scatter(prfROIData.(fnRoi{roi_idx}).x,prfROIData.(fnRoi{roi_idx}).y,[],c,'.');  hold on;
    scatter(prfROIData.(fnRoi{roi_idx}).x_smoothed,prfROIData.(fnRoi{roi_idx}).y_smoothed,[],c_sm,'.');
    
    title(fnRoi{roi_idx})
    axis square
    xlim([-10 10])
    ylim([-10 10])
    
    
    legend({'orig', 'smoothed'},'Location','northeastoutside')
    
end
print(fH8, fullfile(saveDir,'prf_center_distribution'), '-dpng');

%% -----------------------------------------------
% Visualize pRF parameters on brainstorm surface
%-----------------------------------------------
if opt.surfVisualize
    close all;
    
    % Get directory to save images
    saveDir = fullfile(dirPth.fmri.saveDataPth_prfBS, 'figs');
    if ~exist(saveDir, 'dir'); mkdir(saveDir); end
    
    
    mprf_VisualizeDataOnBrainstormSurface(dirPth,surfacesToLoad,saveDir);
    
end


