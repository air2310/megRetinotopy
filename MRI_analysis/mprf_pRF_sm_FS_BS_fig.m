function mprf_pRF_sm_FS_BS_fig(dirPth,opt)
% Function to make a variety of figures after running the function:
%       mprf_pRF_sm_FS_BS()
% Figures include (for different rois and original vs smoothed data)
% 1) distribution of variance explained, eccentricity, sigma, beta and
% recomputed beta values for all Brainstorm vertices
% 2) prf size vs eccentricity
% 3) x vs y position
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
box off; axis square;
set(gca, 'TickDir', 'out', 'FontSize', 14)
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
box off; axis square;
set(gca, 'TickDir', 'out', 'FontSize', 14)
print(fH2, fullfile(saveDir,'voxels_used_for_smoothing'), '-dpng');


%% ROIs

% load wang rois
combineHemi = false;
combineRois = true;
roiData = loadWangROIs(roiBSDir, 'pial.wang2015_atlas', combineHemi, combineRois);

% construct 2 masks: one for Variance Explained above 10%, one with
% eccentricity below 10 deg.
vemask    = (prfData.varexplained > opt.mri.varExplThresh(1) & prfData.varexplained < opt.mri.varExplThresh(2));
eccenmask = (prfData.eccentricity > opt.mri.eccThresh(1) & prfData.eccentricity < opt.mri.eccThresh(2));

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
            %data = data(vemask(roiMask));
            %prfROIData.(fnRoi{r}).(fnPrf{p}) = data;
            
            data = data((vemask(roiMask) & eccenmask(roiMask)));
            
            prfROIData.(fnRoi{r}).(fnPrf{p}) = data;
            
        elseif (strcmp(thisFieldName,'beta') || strcmp(thisFieldName,'recomp_beta'))
            data = data(roiMask);
            thresh = prctile(data, opt.mri.betaPrctileThresh);
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
    box off; axis square;
    set(gca, 'TickDir', 'out', 'FontSize', 14)
end
print(fH3, fullfile(saveDir,'pRF_size_eccentricity'), '-dpng');


% prf size vs ecc for all rois
roiName_allROI = 'allROIs';
fH201 = figure(201); clf; set(gcf, 'Color', 'w', 'Position',[675,384,763,590],'Name', 'pRF size vs ecc all rois');
c_orig 	 = [0.5 1 0];
c_smooth = [0 0.5 1];
scatter(prfROIData.(roiName_allROI).eccentricity, prfROIData.(roiName_allROI).sigma,[],c_orig,'*');  hold on;
scatter(prfROIData.(roiName_allROI).eccentricity_smoothed, prfROIData.(roiName_allROI).sigma_smoothed,[],c_smooth,'*');
title(roiName_allROI, 'FontSize', 15)
legend('original','smoothed','Location','NorthWest')
legend('Location','NorthWest')
ylim([0 opt.mri.eccThresh(2)]); xlim([0 opt.mri.eccThresh(2)]); axis square
ylabel('PRF size (deg)'); xlabel('PRF eccen (deg');
set(gca, 'FontSize', 20, 'TickDir', 'out','LineWidth',3); box off;
print(fH201, fullfile(saveDir,'pRF_size_eccentricity_allROIs'), '-dpng');


% Fit a line to the distribution of points and plot them all in one figure, 
fH31 = figure(31); clf; set(gcf, 'Color', 'w', 'Position', [675,384,763,590], 'Name', 'pRF fit for all rois');

roi_colors = [0.5 0.5 0.5; 1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1; 0.75 0.75 0; 0 0.75 0.75; 0.75 0 0.75...
              ; 0.75 0.75 0.75; 0.75 0 0; 0 0 0.75; 0.75 1 0.75; 1 0 0.75; 0.75 0 1];

roiToPlot = {'V1','V2','V3','IPS','LO','TO','VO','PHC','FEF','SPL1','hV4','allROIs'};          
numRoi    = length(roiToPlot); 
count = 1;
for roiIdx = 1:numRoi
    x    = prfROIData.(roiToPlot{roiIdx}).eccentricity_smoothed;
    y    = prfROIData.(roiToPlot{roiIdx}).sigma_smoothed;
    w    = prfROIData.(roiToPlot{roiIdx}).varexplained;
    xfit = linspace(opt.mri.eccThresh(1),opt.mri.eccThresh(2),20)';
    
    if ~isempty(y)
        [yfit,stats] = mprf_fit(x,y,w,xfit);
        plot(xfit,yfit,'color',roi_colors(roiIdx,:),'LineWidth',4);hold on;
        legend(roiToPlot,'Location','bestoutside')
        xlabel('eccentricity (deg)'); ylabel('pRF size (deg)');
        ylim([0 10]);
        lg{count} = roiToPlot{roiIdx};
        %    xlim(opt.xaxislim);
        set(gca, 'FontSize', 20, 'TickDir','out','LineWidth',3); box off
        count = count+1;
    end
end
legend(lg,'Location','bestoutside');
print(fH31, fullfile(saveDir,'prf_sig_ecc_all_rois'), '-dpng');




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
    box off; axis square;
    set(gca, 'TickDir', 'out', 'FontSize', 14)
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
    box off; axis square;
    set(gca, 'TickDir', 'out', 'FontSize', 14)
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
    box off; axis square;
    set(gca, 'TickDir', 'out', 'FontSize', 14)
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
    box off; axis square;
    set(gca, 'TickDir', 'out', 'FontSize', 14)
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
    set(gca, 'TickDir', 'out', 'FontSize', 14)
    
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


