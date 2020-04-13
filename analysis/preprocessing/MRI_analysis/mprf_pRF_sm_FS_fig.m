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
    tmp = strsplit(curLHFile, '.');
    prfName = tmp{2};
    
    % Find the corresponding rh file:
    curRHFile = ['rh.' prfName];
    
    % load and concatenate:
    bothFSHemiData = [read_curv(fullfile(prfFSDir,curLHFile));...
        read_curv(fullfile(prfFSDir,curRHFile))];
    
    prfDataFS.(prfName) = bothFSHemiData;
    
end

%% Load FS ROIs

if opt.roi.roimrvToFS
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
  
        roiDataFS.(roiName) = find(bothHemiROI);
    end
 
else
    combineHemi = true;
    combineRois = true;
    roiDataFS = loadWangROIs(dirPth.fs.surfPth, '*.wang2015_atlas.mgz', combineHemi, combineRois);
end

% Sample prf data for different visual areas
vemask     = (prfDataFS.varexplained > opt.mri.varExplThresh(1) & prfDataFS.varexplained < opt.mri.varExplThresh(2));
eccenmask  = (prfDataFS.eccentricity > opt.mri.eccThresh(1) & prfDataFS.eccentricity < opt.mri.eccThresh(2));
prfROIData = struct();
fnData     = fieldnames(prfDataFS);
fnRoi      = fieldnames(roiDataFS);
numROIs = length(fnRoi);

for p = 1:length(fnData)

    for r = 1:numROIs

        % Get prf data
        thisFieldName = fnData{p};
        data = prfDataFS.(fnData{p});
        
        % Get roi mask
        roimask = roiDataFS.(fnRoi{r})>0;
        
        % Combine masks
        allMasks = (vemask&eccenmask&roimask);
        
        % Mask data
        prfROIData.(fnRoi{r}).(fnData{p}) = data(allMasks);
        
    end
        
end



% Subplot variables
nRow = 3;
nCol = ceil(numROIs/nRow);
scr_sz = get(0, 'screensize');

% Variance explained (all voxels, whole brain)
%-------------------------------
scr_sz = get(0, 'screensize');
fH1 = figure(101); clf; set(gcf, 'Color', 'w', 'Position', [1 1 scr_sz(3)/2, scr_sz(4)/2]);
h1 = histogram(prfDataFS.varexplained,100); hold on;
h2 = histogram(prfDataFS.varexplained(vemask),100);
h1.FaceColor = 'g'; h1.EdgeColor = 'w';
h2.FaceColor = 'r'; h2.EdgeColor = 'w';
title('Variance explained of pRF model fits on FS surface (all vertices, or (ve>0.1)')
xlabel('Variance explained'); ylabel('Frequency')
set(gca, 'FontSize', 15, 'TickDir', 'out'); box off;
print(fH1, fullfile(saveDir,'variance_explained_selection'), '-dpng');

% pRF size vs eccentricity
%----------------------
fH3 = figure(103); clf; set(gcf, 'Color', 'w', 'Position', scr_sz);
c = [0.5 1 0]; 
c_sm = [0 0.5 1];
title('pRF size vs ecc')
for roiIdx = 1:numROIs
    subplot(nRow,nCol,roiIdx);
    
    scatter(prfROIData.(fnRoi{roiIdx}).eccentricity, prfROIData.(fnRoi{roiIdx}).sigma,[],c,'*');  hold on;
    scatter(prfROIData.(fnRoi{roiIdx}).eccentricity_smoothed,prfROIData.(fnRoi{roiIdx}).sigma_smoothed,[],c_sm,'*');
    
    title(fnRoi{roiIdx}, 'FontSize', 15)
    legend('Original','Smoothed','Location','NorthWest')
    ylim([0 opt.mri.eccThresh(2)]); xlim([0 opt.mri.eccThresh(2)]); axis square
    ylabel('PRF sigma (deg)'); xlabel('PRF eccen (deg)');
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH3, fullfile(saveDir,'pRF_size_eccentricity'), '-dpng');

% PRF size vs ecc across all rois
%----------------------
roiName_allROI = fnRoi{1};
fH4 = figure(104); clf; set(gcf, 'Color', 'w', 'Position',[675,384,763,590],'Name', 'pRF size vs ecc all rois');
scatter(prfROIData.(roiName_allROI).eccentricity, prfROIData.(roiName_allROI).sigma,[],c,'*');  hold on;
scatter(prfROIData.(roiName_allROI).eccentricity_smoothed, prfROIData.(roiName_allROI).sigma_smoothed,[],c_sm,'*');

title('allROIs', 'FontSize', 15)
legend('Original','Smoothed','Location','NorthWest'); legend boxoff;
ylim([0 opt.mri.eccThresh(2)]); xlim([0 opt.mri.eccThresh(2)]); axis square
ylabel('PRF size (deg)'); xlabel('PRF eccen (deg');
set(gca, 'FontSize', 20, 'TickDir', 'out','LineWidth',2); box off;
print(fH4, fullfile(saveDir,'pRF_size_eccentricity_allROIs'), '-dpng');

%% Plot line fits with size vs eccen data
scr_sz = get(0,'screensize');
fH5 = figure(105); clf; set(gcf, 'Color', 'w', 'Position', scr_sz, 'Name', 'pRF size vs ecc with fits');

roisToPlot = cellfind(fnRoi, 'V1'):numROIs;
roi_colors =  hsv(length(roisToPlot));

for roiIdx = 1:length(roisToPlot)
    subplot(nRow,nCol,roiIdx);

    if ~isempty(prfROIData.(fnRoi{roisToPlot(roiIdx)}).varexplained)
        x    = prfROIData.(fnRoi{roisToPlot(roiIdx)}).eccentricity_smoothed;
        y    = prfROIData.(fnRoi{roisToPlot(roiIdx)}).sigma_smoothed;
        w    = prfROIData.(fnRoi{roisToPlot(roiIdx)}).varexplained;
        xfit = linspace(opt.mri.eccThresh(1),opt.mri.eccThresh(2),20)';
        yfit = mprf_fit(x,y,w,xfit);

        scatter(x, y, [], roi_colors(roiIdx,:),'*'); hold on;
        plot(xfit,yfit,'color','k','LineWidth',4); 
    end
    title(fnRoi{roisToPlot(roiIdx)}, 'FontSize', 15)
    ylim([0 opt.mri.eccThresh(2)]); xlim([0 opt.mri.eccThresh(2)]); axis square
    ylabel('PRF size (deg)'); xlabel('PRF eccen (deg');
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH5, fullfile(saveDir,'pRF_size_eccentricity_smoothed_roi_wFits'), '-dpng');


%% Fit a line to the distribution of points and plot them all in one figure, 
fH6 = figure(106); clf; set(gcf, 'Color', 'w', 'Position', [1 1 scr_sz(3)/2 scr_sz(4)/2], 'Name', 'pRF fit for all rois');

for roiIdx = 1:length(roisToPlot)
    if ~isempty(prfROIData.(fnRoi{roisToPlot(roiIdx)}).varexplained)
        x    = prfROIData.(fnRoi{roisToPlot(roiIdx)}).eccentricity_smoothed;
        y    = prfROIData.(fnRoi{roisToPlot(roiIdx)}).sigma_smoothed;
        w    = prfROIData.(fnRoi{roisToPlot(roiIdx)}).varexplained;

        xfit = linspace(opt.mri.eccThresh(1),opt.mri.eccThresh(2),20)';
        yfit = mprf_fit(x,y,w,xfit);

        plot(xfit,yfit,'color',roi_colors(roiIdx,:),'LineWidth',4);hold on;
    end 
end
legend(fnRoi{roisToPlot},'location','bestoutside'); title('pRF size vs eccentricity (linear fit)');
legend boxoff;
ylim([0 opt.mri.eccThresh(2)]); xlim([0 opt.mri.eccThresh(2)]); axis square
ylabel('PRF size (deg)'); xlabel('PRF eccen (deg');
set(gca, 'FontSize', 14, 'TickDir', 'out'); box off
print(fH6, fullfile(saveDir,'pRF_size_eccentricity_smoothed_fits_all_rois'), '-dpng');


%% Histograms

% Sigma histogram
%-----------------------
fH7 = figure(107); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);

for roiIdx = 1:numROIs
    subplot(nRow,nCol,roiIdx);
    hist(prfROIData.(fnRoi{roiIdx}).sigma,100)
    
    title(fnRoi{roiIdx}, 'FontSize', 15)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    xlabel('Sigma'); ylabel('Frequency')
    box off;
end
print(fH7, fullfile(saveDir,'sigma'), '-dpng');

% Smoothed sigma histogram
%-----------------------
fH8 = figure(108); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);

for roiIdx = 1:numROIs
    subplot(nRow,nCol,roiIdx);
    hist(prfROIData.(fnRoi{roiIdx}).sigma_smoothed,100)
    
    title(fnRoi{roiIdx}, 'FontSize', 15)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    xlabel('Sigma smoothed'); ylabel('Frequency')
    box off;
end
print(fH5, fullfile(saveDir,'sigma_smoothed'), '-dpng');

% Beta histogram
%-----------------------
fH9 = figure(109); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);

for roiIdx = 1:numROIs
    subplot(nRow,nCol,roiIdx);
    hist(prfROIData.(fnRoi{roiIdx}).beta,100)
    
    title(fnRoi{roiIdx}, 'FontSize', 15)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    box off;
    xlabel('Beta'); ylabel('Frequency')
end
print(fH9, fullfile(saveDir,'beta'), '-dpng');

% recomputed beta histogram
%-----------------------
fH10 = figure(110); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
for roiIdx = 1:numROIs
    subplot(nRow,nCol,roiIdx);
    hist(prfROIData.(fnRoi{roiIdx}).recomp_beta,100)
    title(fnRoi{roiIdx}, 'FontSize', 15)
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    box off;
    xlabel('Recomp Beta'); ylabel('Frequency')
end
print(fH10, fullfile(saveDir,'recomp_beta'), '-dpng');

fH11 = figure(111); set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080]);
c = [0.5 1 0];
c_sm = [0 0.5 1];
title('pRF center distribution')
for roiIdx = 1:numROIs
    subplot(nRow,nCol,roiIdx);
    scatter(prfROIData.(fnRoi{roiIdx}).x,prfROIData.(fnRoi{roiIdx}).y,[],c,'.');  hold on;
    scatter(prfROIData.(fnRoi{roiIdx}).x_smoothed,prfROIData.(fnRoi{roiIdx}).y_smoothed,[],c_sm,'.');
    %axis image
    xlim([-10 10])
    ylim([-10 10])
    set(gca, 'TickDir', 'out','FontSize',15);
    box off;
    title(fnRoi{roiIdx}, 'FontSize', 15)
    legend({'Orig', 'Smoothed'}, 'Location','NorthWest')
    
    
end
print(fH11, fullfile(saveDir,'prf_center_distribution'), '-dpng');

%% --------------------------------------------------------------------
%   Load mrMesh and display Wang Rois and prf parameters on the mesh
%  --------------------------------------------------------------------

if opt.surfVisualize
    close all;
    
    % Get directory to save images
    saveDir = fullfile(dirPth.fmri.saveDataPth_prfFS, 'figs');
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
            
            mprf_VisualizeDataOnFreesurferSurface(dirPth, cur_surf_to_load, saveDir, opt);
            
            mprf_VisualizeRoiOnFreesurferSurface(dirPth, cur_surf_to_load,saveDir, opt);
        end
    end
end

end