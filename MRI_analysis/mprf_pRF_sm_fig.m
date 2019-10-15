function mprf_pRF_sm_fig(dirPth, opt)
% Function to visualize prf parameters as histograms and on mrVista surface
% after completing the function mprf_pRF_sm().
%   1) distribution of variance explained, beta, prf size, prf position
%       values for all voxels in visual ROIs and across all ROIs
% 	2) prf size vs eccentricity (different rois)
%                               (original and smoothed)
%   3) mrVista surface plots (original and smoothed)
%
% INPUTS:
%   subjID      :   subject name (string)
%   dirPth      :   paths locating subject's data and files (struct, see loadPaths.m)
%
% OUTPUTS:
%   none
%
%see also mprf_pRF_sm

%% ------------
%  0. Get file paths and load mrVista retinotopy data
% -------------
anatDir         = fullfile(dirPth.fmri.mrvPth, '3DAnatomy');
grayCoordsDir   = fullfile(dirPth.fmri.mrvPth, 'Gray');
roiDir          = fullfile(anatDir, 'ROIs');

% Get original prf parameters from mat-file
load(fullfile(dirPth.fmri.saveDataPth_prfMrv, 'mat', 'exported_prf_params.mat'), 'prf_par_exp');
prf = prf_par_exp; clear prf_par_exp;

load(fullfile(grayCoordsDir,'coords.mat'), 'coords');

% Get directory to save images
saveDir = fullfile(dirPth.fmri.saveDataPth_prfMrv, 'figs');
if ~exist(saveDir, 'dir'); mkdir(saveDir); end

%% ------------------------------------------------
%   Separate prf params for different ROIs
%  ------------------------------------------------

% Load Wang Atlas ROIs compatible with mrVista, combine
% dorsal/ventral/etc rois
roisToCombine = {'V1','V2','V3',...
    'IPS','LO','TO','VO', 'PHC'};

d = dir(fullfile(roiDir, 'WangAtlas_*'));

% Get rid of combined files
toRemove = zeros(1,length(d));
for nn = 1:length(d)
    
    curMatFile = d(nn).name;
    result = regexp(curMatFile, 'Combined', 'match');
    if ~isempty(result)
        toRemove(nn) = 1;
    end
end
d(logical(toRemove)) = [];

% Put struct in a table
t = struct2table(d);

clear roiName roiLoc

for rc = 1:length(roisToCombine)
    match = contains(t.name,roisToCombine(rc));
    
    matchIdx = find(match);
    tmpData = [];
    for ii = 1:length(matchIdx)
        load(fullfile(roiDir, t.name{matchIdx(ii)}));
        [~, indices] = intersect(coords', ROI.coords', 'rows' );
        tmpData = [tmpData; indices];
    end
    
    roiName(rc) = roisToCombine(rc);
    roiLoc{rc} = unique(tmpData, 'rows');
    
end

singleRois = ~contains(t.name,roisToCombine);
for sr = find(singleRois)'
    
    % Remove irrelevant parts of filename
    tmp = strsplit(t.name{sr}, {'.','_'});
    
    % Add ROI name and location to list
    roiName(end+1) = tmp(2);
    
    load(fullfile(roiDir, t.name{sr}));
    [~, indices] = intersect(coords', ROI.coords', 'rows');
    roiLoc{end+1} = indices;
end

% Apply roi masks to pRF parameters
numRoi = length(roiName);
prfROIData = struct();
fn = fieldnames(prf);

vemask    = (prf.varexplained > opt.mri.varExplThresh(1) & prf.varexplained < opt.mri.varExplThresh(2));
eccenmask = (prf.eccentricity > opt.mri.eccThresh(1) & prf.eccentricity < opt.mri.eccThresh(2));

% percentage of voxels removed overall out of all the voxels (might be
% different for 004 and 040 because of the small FOV in these two subjects)
allRoi = 0;             
mask_roi_allRoi = 0; 
mask_all      = sum(vemask & eccenmask);
excluded.all  = (1-mask_all/size(prf.varexplained,2)).*100;
fprintf('\n percentage of total voxels excluded = %.2f \n',excluded.all);
for roiIdx = 1:numRoi 
    curRoi                  = roiName{roiIdx};
    mask_roi_curRoi         = (vemask(roiLoc{roiIdx}) & eccenmask(roiLoc{roiIdx}));
    
    allRoi                  = allRoi + size(roiLoc{roiIdx},1);
    mask_roi_allRoi         = mask_roi_allRoi + sum(mask_roi_curRoi);
    
    excluded.(curRoi)       = (1-sum(mask_roi_curRoi)/size(roiLoc{roiIdx},1)).*100;
    excluded_rois(:,roiIdx)  = excluded.(curRoi);  
        
    fprintf('\n percentage of voxels excluded from roi %s = %.2f \n',curRoi,excluded.(curRoi)); 
end

excluded.allRoi = (1-mask_roi_allRoi/allRoi).*100;
excluded_rois(:,end+1)  = excluded.allRoi;


for p = 1:length(fn)
    
    allroitmp = [];
    % determine the percentage of voxels removed after the ve and
    % eccentricity threshold
   
    
    for roiIdx = 1:numRoi
        
        thisFieldName = fn{p};
        data = prf.(fn{p});
        
        if strcmp(thisFieldName,'varexplained')
            data = data(roiLoc{roiIdx});
            %prfROIData.(roiName{roiIdx}).(fn{p}) = data(data >
            %opt.mri.varExplThresh(1) & data < opt.mri.varExplThresh(2));
            
            data = data((vemask(roiLoc{roiIdx}) & eccenmask(roiLoc{roiIdx})));
            
            prfROIData.(roiName{roiIdx}).(fn{p}) = data;
            
        elseif (strcmp(thisFieldName,'beta') || strcmp(thisFieldName,'recomp_beta'))
            data = data(roiLoc{roiIdx});
            thresh = prctile(data, opt.mri.betaPrctileThresh);
            betamask = ((data > thresh(1)) & (data < thresh(2)));
            prfROIData.(roiName{roiIdx}).(fn{p}) = data(betamask);
        else
            
            % Mask
            data = data(roiLoc{roiIdx});
            data = data((vemask(roiLoc{roiIdx}) & eccenmask(roiLoc{roiIdx})));
            
            prfROIData.(roiName{roiIdx}).(fn{p}) = data;
        end
        allroitmp  = [allroitmp, prfROIData.(roiName{roiIdx}).(fn{p})];
    end
    prfROIData.allROIs.(thisFieldName) = allroitmp;
end

roiName(end+1) = {'allROIs'};
numRois = length(roiName);


% make a bar graph of the voxels excluded
fH100 = figure(100); clf; set(gcf, 'Color', 'w', 'Position', [41,1,1876,973]);
bar(excluded_rois);
title('Voxels excluded')
xlabel('visual areas'); ylabel('voxels excluded (%)')
set(gca, 'FontSize', 30, 'TickDir', 'out','xTickLabel',roiName); box off;
ylim([0 100]);
print(fH100, fullfile(saveDir,'variance_explained'), '-dpng');



%% ----------------------
%  1. Visualize distribution of prf parameters for all voxels (whole brain)
%  ----------------------

% Variance explained (all voxels, whole brain)
fH1 = figure(1); clf; set(gcf, 'Color', 'w', 'Position', [41,1,1876,973]);
hist(prf.varexplained,100);
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
title('Variance explained of pRF model fits (all vertices, whole brain)')
xlabel('Variance explained'); ylabel('Frequency')
set(gca, 'FontSize', 15, 'TickDir', 'out'); box off;
print(fH1, fullfile(saveDir,'variance_explained'), '-dpng');

% median variance explained for each roi
med_ve = median(prfROIData.allROIs.varexplained);
fprintf('median variance explained = %.2f',med_ve);

%% ------------------------------------------------
%  2. Visualize prf params for different ROIs
%  ------------------------------------------------

nrows = 2;
ncols = ceil(numRois/nrows);
nbins = 50;

% Plot pRF size vs ecc
%-----------------------
fH2 = figure(2); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'pRF size vs ecc');
c_orig 	 = [0.5 1 0];
c_smooth = [0 0.5 1];

for roiIdx = 1:numRois
    subplot(nrows,ncols,roiIdx);
    scatter(prfROIData.(roiName{roiIdx}).eccentricity, prfROIData.(roiName{roiIdx}).sigma,[],c_orig,'*');  hold on;
    scatter(prfROIData.(roiName{roiIdx}).eccentricity_smoothed, prfROIData.(roiName{roiIdx}).sigma_smoothed,[],c_smooth,'*');
    
    title(roiName{roiIdx}, 'FontSize', 15)
    legend('original','smoothed','Location','NorthWest')
    legend('Location','NorthWest')
    ylim([0 opt.mri.eccThresh(2)]); xlim([0 opt.mri.eccThresh(2)]); axis square
    ylabel('PRF size (deg)'); xlabel('PRF eccen (deg');
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH2, fullfile(saveDir,'pRF_size_eccentricity'), '-dpng');

% prf size vs ecc for all rois
roiName_allROI = roiName{end};
fH21 = figure(21); clf; set(gcf, 'Color', 'w', 'Position',[675,384,763,590],'Name', 'pRF size vs ecc all rois');
c_orig 	 = [0.5 1 0];
c_smooth = [0 0.5 1];
scatter(prfROIData.(roiName_allROI).eccentricity, prfROIData.(roiName{roiIdx}).sigma,[],c_orig,'*');  hold on;
scatter(prfROIData.(roiName_allROI).eccentricity_smoothed, prfROIData.(roiName{roiIdx}).sigma_smoothed,[],c_smooth,'*');
title(roiName_allROI, 'FontSize', 15)
legend('original','smoothed','Location','NorthWest')
legend('Location','NorthWest')
ylim([0 opt.mri.eccThresh(2)]); xlim([0 opt.mri.eccThresh(2)]); axis square
ylabel('PRF size (deg)'); xlabel('PRF eccen (deg');
set(gca, 'FontSize', 20, 'TickDir', 'out','LineWidth',3); box off;
print(fH21, fullfile(saveDir,'pRF_size_eccentricity_allROIs'), '-dpng');



fH3 = figure(3); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'Variance explained');
% Variance explained of the model fit split into rois
for roiIdx = 1:numRois
    subplot(nrows,ncols,roiIdx);
    scatter(prfROIData.(roiName{roiIdx}).eccentricity, prfROIData.(roiName{roiIdx}).varexplained,[],c_orig,'*');  hold on;
    
    title(roiName{roiIdx}, 'FontSize', 15)
    legend('original','smoothed','Location','NorthWest')
    legend('Location','NorthWest')
    ylim([0 1]); xlim([0 opt.mri.eccThresh(2)]); axis square
    ylabel('variance explained (deg)'); xlabel('PRF eccen (deg');
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH3, fullfile(saveDir,'variance_explained'), '-dpng');


% Fit a line to the distribution of points and plot them all in one figure, 
fH4 = figure(4); clf; set(gcf, 'Color', 'w', 'Position', [675,384,763,590], 'Name', 'pRF fit for all rois');

roi_colors = [0.5 0.5 0.5; 1 0.5 0.5; 0.5 1 0.5; 0.5 0.5 1; 0.75 0.75 0; 0 0.75 0.75; 0.75 0 0.75...
              ; 0.75 0.75 0.75; 0.75 0 0; 0 0 0.75; 0.75 1 0.75; 1 0 0.75; 0.75 0 1];
count = 1;

for roiIdx = 1:numRois
    x    = prfROIData.(roiName{roiIdx}).eccentricity_smoothed;
    y    = prfROIData.(roiName{roiIdx}).sigma_smoothed;
    w    = prfROIData.(roiName{roiIdx}).varexplained;
    xfit = linspace(opt.mri.eccThresh(1),opt.mri.eccThresh(2),20)';
    
    if ~isempty(y)
        [yfit,stats] = mprf_fit(x,y,w,xfit);
        plot(xfit,yfit,'color',roi_colors(roiIdx,:),'LineWidth',4);hold on;
        
        lg{count} = roiName{roiIdx};
        
        xlabel('eccentricity (deg)'); ylabel('pRF size (deg)');
        %    ylim(opt.yaxislim);
        %    xlim(opt.xaxislim);
        set(gca, 'FontSize', 20, 'TickDir','out','LineWidth',3); box off
        hold on;
        count = count+1;
    end

end
legend(lg,'Location','bestoutside');
print(fH4, fullfile(saveDir,'prf_sig_ecc_all_rois'), '-dpng');



% Sigma histogram
%-----------------------
fH5 = figure(5); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'Sigma');
for roiIdx = 1:numRois
    subplot(nrows,ncols,roiIdx);
    hist(prfROIData.(roiName{roiIdx}).sigma,nbins)
    
    title(roiName{roiIdx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    xlim([0 opt.mri.eccThresh(2)]); box off;
    ylabel('Frequency'); xlabel('PRF size (deg)')
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH5, fullfile(saveDir,'sigma'), '-dpng');

% Smoothed sigma histogram
%-----------------------
fH6 = figure(6); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name','Sigma smoothed');
for roiIdx = 1:numRois
    subplot(nrows,ncols,roiIdx);
    
    hist(prfROIData.(roiName{roiIdx}).sigma_smoothed,nbins)
    title(roiName{roiIdx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    xlim([0 opt.mri.eccThresh(2)]); box off;
    ylabel('Frequency'); xlabel('PRF smoothed size (deg)')
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH6, fullfile(saveDir,'sigma_smoothed'), '-dpng');

% Beta histogram
%-----------------------
fH7 = figure(7); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name','Beta')
maxBeta = prctile(prfROIData.allROIs.beta,[opt.mri.betaPrctileThresh(1) opt.mri.betaPrctileThresh(2)]);

for roiIdx = 1:numRois
    subplot(nrows,ncols,roiIdx);
    dataToPlot = prfROIData.(roiName{roiIdx}).beta;
    betaMask = dataToPlot<maxBeta(2);
    hist(dataToPlot(betaMask),nbins)
    title(roiName{roiIdx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    box off;
    ylabel('Frequency'); xlabel('PRF beta (a.u.)')
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH7, fullfile(saveDir,'beta'), '-dpng');

% recomputed beta histogram
%-----------------------
fH8 = figure(8); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'Recomputed beta')
maxBeta = prctile(prfROIData.allROIs.recomp_beta,[opt.mri.betaPrctileThresh(1) opt.mri.betaPrctileThresh(2)]);
for roiIdx = 1:numRois
    subplot(nrows,ncols,roiIdx);
    dataToPlot = prfROIData.(roiName{roiIdx}).recomp_beta;
    betaMask = dataToPlot<maxBeta(2);
    hist(dataToPlot(betaMask),nbins)
    title(roiName{roiIdx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    box off;
    ylabel('Frequency'); xlabel('PRF Recomp beta (a.u.)')
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH8, fullfile(saveDir,'recomp_beta'), '-dpng');


fH9 = figure(9); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920/2   1080], 'Name', 'pRF center distribution');
% polar plot params
params.grid          = 'on';
params.line          = 'off';
params.gridColor     = [0.7,0.7,0.7];
params.fontSize      = 20;
params.symbol        = 'o';
params.size          = 1;
params.color         = 'w';
params.fillColor     = 'w';
params.maxAmp        = 15;
params.ringTicks     = [0:5:10];
params.gridLineWidth = 0.01;
params.lineWidth     = 0.2;

for roiIdx = 1:numRois
    
    subplot(nrows,ncols,roiIdx);
    title(roiName{roiIdx})
    polarPlot(0,params); hold on;
    
    x = prfROIData.(roiName{roiIdx}).x;
    y = prfROIData.(roiName{roiIdx}).y;
    x_smoothed = prfROIData.(roiName{roiIdx}).x_smoothed;
    y_smoothed  = prfROIData.(roiName{roiIdx}).y_smoothed;
    
    plot(x,y,'o','MarkerFaceColor',c_orig,'MarkerSize',2); hold all;
    plot(x_smoothed,y_smoothed,'o','MarkerFaceColor',c_smooth,'MarkerSize',2);
    
end
print(fH9, fullfile(saveDir,'prf_center_distribution'), '-dpng');



    
%% --------------------------------------------------------------------
%   Load mrMesh and display Wang Rois and prf parameters on the mesh
%  --------------------------------------------------------------------

if opt.surfVisualize
    
    % Go to vista session and open a mrVista Gray window
    cd(dirPth.fmri.mrvPth)
    vw = mrVista('3');
    
    % Load rh and lh mesh
    mesh1 = fullfile(anatDir, 'Left', '3DMeshes', 'Left_inflated.mat');
    mesh2 = fullfile(anatDir, 'Right', '3DMeshes', 'Right_inflated.mat');
    [vw, OK] = meshLoad(vw, mesh1, 1); if ~OK, error('Mesh server failure'); end
    [vw, OK] = meshLoad(vw, mesh2, 1); if ~OK, error('Mesh server failure'); end
    
    % Get prf parameters saved as nifti's
    prfParams = {'eccentricity', 'eccentricity_smoothed', 'polar_angle', 'polar_angle_smoothed', ...
        'sigma', 'sigma_smoothed', 'varexplained', 'beta', 'recomp_beta'};
    
    % Hide rois in gray view when loading
    vw = viewSet(vw, 'hide gray rois', true);
    
    % Load and draw Wang ROIs
    for idx = 1:length(t.name)
        roiFile = sprintf('%s',t.name{idx});
        vw = loadROI(vw, roiFile);
        fprintf('(%s): Loaded ROI: %s \n', mfilename, roiFile)
    end
    
    vw = viewSet(vw, 'ROI draw method', 'perimeter');
    vw = refreshScreen(vw);
    vw = meshUpdateAll(vw);
    
    % Get list of viewpoints and meshes when saving different images
    viewList={'back','left','right','bottom','top'};
    viewVectors={[pi -pi/2 0],[pi 0 0],[0 0 pi],[pi/2 -pi/2 0],[-pi/2 -pi/2 0]};
    mshToPrint = viewGet(vw,'allmeshes');
    idList = [mshToPrint{1}.id, mshToPrint{2}.id];
    idName = {mshToPrint{1}.name, mshToPrint{2}.name};
    
    % Load and draw prf parameters
    fH = figure('Color', 'w');
    
    for ii = 1:length(prfParams)
        
        fprintf('(%s):  Visualizing %s on mrVista surface \n', mfilename, prfParams{ii})
        
        prfParamNifti =  fullfile(dirPth.fmri.saveDataPth_prfMrv,'nifti', sprintf('%s.nii.gz', prfParams{ii}));
        
        data_type = 'Averages';
        vw = viewSet(vw,'curdt',data_type);  
        vw = viewSet(vw,'displaymode','map');
        vw = loadParameterMap(vw,prfParamNifti);
        vw = viewSet(vw,'cothresh',opt.mri.varExplThresh(1));
        
        
        % Add smoothing/ve mask?
        %smoothing_mask_tmp = (vemask);
        %smoothing_mask     = nan(size(smoothing_mask_tmp));
        %smoothing_mask(smoothing_mask_tmp==1) = 1;
        vw.map{1} = vw.map{1}; %.*smoothing_mask;
        %vw.map{1} = vw.map{1};
        
        switch prfParams{ii}
            case {'polar_angle', 'polar_angle_smoothed'}
                cAxisLim = [eps 2*pi];
                vw = viewSet(vw, 'mapwin', cAxisLim);
                vw = viewSet(vw, 'mapclip', cAxisLim);
                vw.ui.mapMode = setColormap(vw.ui.mapMode, 'hsvCmap');
                
            case {'eccentricity', 'eccentricity_smoothed'}
                % Set colormap and limits
                cAxisLim = [eps opt.mri.eccThresh(2)];
                vw.ui.mapMode = setColormap(vw.ui.mapMode, 'hsvTbCmap');
                vw = viewSet(vw, 'mapwin', cAxisLim);
                vw = viewSet(vw, 'mapclip', cAxisLim);
                
            case {'beta', 'recomp_beta'}
                if opt.mri.betaPrctileThresh
                    maxBetaCmap = prctile(vw.map{1},[opt.mri.betaPrctileThresh(1) opt.mri.betaPrctileThresh(2)]);
                    cAxisLim = [eps maxBetaCmap];
                    vw = viewSet(vw, 'mapwin', cAxisLim);
                    vw = viewSet(vw, 'mapclip', cAxisLim);
                end
        end
        
        % Update views
        vw = refreshScreen(vw);
        vw = meshUpdateAll(vw);
        
        % Copy the mesh to a Matlab figure
        for thisID = 1:length(idList)
            m=mrmSet(vw.mesh{thisID}, 'actor');
            
            for thisView=1:length(viewList)
                cam.actor=0;
                cam.rotation = rotationMatrix3d(viewVectors{thisView});
                mrMesh('localhost',m.id,'set',cam);
                
                figure(fH); clf;
                imagesc(mrmGet(m, 'screenshot')/255); caxis(cAxisLim); colorbar; axis image; axis off;
                if opt.saveFig print(fH, fullfile(saveDir,sprintf('%s_%s_%s',idName{thisID}, prfParams{ii},viewList{thisView})), '-dpng'); end
                
            end
        end
    end
end


end

