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

%% Print percentage of voxels removed overall out of all the voxels
% (might be different for subjects 004 and 040 because of the small FOV)
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

    fprintf('(%s): percentage of voxels excluded from roi %s = %.2f\n',mfilename,curRoi,excluded.(curRoi));
end

excluded.allRoi = (1-mask_roi_allRoi/allRoi).*100;
excluded_rois(:,end+1)  = excluded.allRoi;

% make a bar graph of the voxels excluded
fH100 = figure(100); clf; set(gcf, 'Color', 'w', 'Position', [41,1,1876,973]);
bar(excluded_rois);
title('Voxels excluded')
xlabel('visual areas'); ylabel('voxels excluded (%)')
set(gca, 'FontSize', 30, 'TickDir', 'out','xTickLabel',roiName); box off;
ylim([0 100]);
print(fH100, fullfile(saveDir,'variance_explained'), '-dpng');


%% Get prf data per ROI
for p = 1:length(fn)
    
    allroitmp = [];   
    
    for roiIdx = 1:numRoi
        % Get prf param name
        thisFieldName = fn{p};
        
        % Select prf data
        data = prf.(fn{p});
        
        % Mask prf data for ROI
        data = data(roiLoc{roiIdx});
            
        data = data((vemask(roiLoc{roiIdx}) & eccenmask(roiLoc{roiIdx})));
        
        % Enter data in struct
        prfROIData.(roiName{roiIdx}).(fn{p}) = data;
                
        % Accumulate across ROIs
        allroitmp  = [allroitmp, prfROIData.(roiName{roiIdx}).(fn{p})];
    end
    prfROIData.allROIs.(thisFieldName) = allroitmp;
end

roiName(end+1) = {'allROIs'};
numRois = length(roiName);
scr_sz = get(0, 'screensize');


%% ----------------------
%  1. Visualize distribution of prf parameters for all voxels (whole brain)
%  ----------------------

% Variance explained (all voxels, whole brain)
fH1 = figure(1); clf; set(gcf, 'Color', 'w', 'Position', [1 1 scr_sz(3)/2, scr_sz(4)/2]);
h1 = histogram(prf.varexplained,100); hold on;
h2 = histogram(prf.varexplained(vemask),100);
h1.FaceColor = 'g'; h1.EdgeColor = 'w';
h2.FaceColor = 'r'; h2.EdgeColor = 'w';
title('Variance explained of pRF model fits (all vertices, or (ve>0.1)')
xlabel('Variance explained'); ylabel('Frequency')
set(gca, 'FontSize', 15, 'TickDir', 'out'); box off;
print(fH1, fullfile(saveDir,'variance_explained'), '-dpng');

% median variance explained for each roi
med_ve = median(prfROIData.allROIs.varexplained);
fprintf('(%s): Median variance explained = %.2f\n',mfilename,med_ve);

%% ------------------------------------------------
%  2. Visualize prf params for different ROIs
%  ------------------------------------------------

nrows = 2;
ncols = ceil(numRois/nrows);
nbins = 50;


roi_colors = hsv(numRois);

% Plot pRF size vs ecc
%-----------------------
fH2 = figure(2); clf; set(gcf, 'Color', 'w', 'Position', scr_sz, 'Name', 'pRF size vs ecc');
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

% prf size vs ecc across all rois
roiName_allROI = roiName{end};
fH3 = figure(3); clf; set(gcf, 'Color', 'w', 'Position',[675,384,763,590],'Name', 'pRF size vs ecc all rois');
scatter(prfROIData.(roiName_allROI).eccentricity, prfROIData.(roiName{roiIdx}).sigma,[],c_orig,'*');  hold on;
scatter(prfROIData.(roiName_allROI).eccentricity_smoothed, prfROIData.(roiName{roiIdx}).sigma_smoothed,[],c_smooth,'*');
title(roiName_allROI, 'FontSize', 15)
legend('original','smoothed','Location','NorthWest'); legend boxoff
ylim([0 opt.mri.eccThresh(2)]); xlim([0 opt.mri.eccThresh(2)]); axis square
ylabel('PRF size (deg)'); xlabel('PRF eccen (deg');
set(gca, 'FontSize', 20, 'TickDir', 'out','LineWidth',2); box off;
print(fH3, fullfile(saveDir,'pRF_size_eccentricity_allROIs'), '-dpng');

%% Variance explained of the model fit split into rois

fH4 = figure(4); clf; set(gcf, 'Color', 'w', 'Position', scr_sz, 'Name', 'Variance explained');
for roiIdx = 1:numRois
    subplot(nrows,ncols,roiIdx);
    scatter(prfROIData.(roiName{roiIdx}).eccentricity, prfROIData.(roiName{roiIdx}).varexplained,[],roi_colors(roiIdx,:),'*');  hold on;
    scatter(prfROIData.(roiName{roiIdx}).eccentricity_smoothed, prfROIData.(roiName{roiIdx}).varexplained,[],roi_colors(roiIdx,:),'*');  hold on;
    
    title(roiName{roiIdx}, 'FontSize', 15)
    legend('original','smoothed','Location','NorthWest')
    legend('Location','NorthWest')
    ylim([0 1]); xlim([0 opt.mri.eccThresh(2)]); axis square
    ylabel('variance explained (deg)'); xlabel('PRF eccen (deg');
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH4, fullfile(saveDir,'variance_explained_selection'), '-dpng');


%% Fit a line to size vs eccen distribution and plot the line fits 
fH4 = figure(4); clf; set(gcf, 'Color', 'w', 'Position', [1 1 scr_sz(3)/2, scr_sz(4)/2], ...
                'Name', 'pRF fit for all rois');

for roiIdx = 1:numRois-1
    x    = prfROIData.(roiName{roiIdx}).eccentricity_smoothed;
    y    = prfROIData.(roiName{roiIdx}).sigma_smoothed;
    w    = prfROIData.(roiName{roiIdx}).varexplained;
    xfit = linspace(opt.mri.eccThresh(1),opt.mri.eccThresh(2),20)';
    
    yfit = mprf_fit(x,y,w,xfit);
    plot(xfit,yfit,'color',roi_colors(roiIdx,:),'LineWidth',2); hold on;
    legend(roiName,'Location','bestoutside'); legend boxoff;
    xlabel('eccentricity (deg)'); ylabel('pRF size (deg)');
    set(gca, 'FontSize', 14, 'TickDir','out','LineWidth',2); box off
    hold on; axis square; title('pRF size vs eccentricity (linear fit)') 

end
print(fH4, fullfile(saveDir,'prf_sig_ecc_all_rois'), '-dpng');

%% Plot line fits with size vs eccen data 
fH5 = figure(5); clf; set(gcf, 'Color', 'w', 'Position', scr_sz, 'Name', 'pRF size vs ecc with fits');

for roiIdx = 1:numRois-1
    subplot(nrows,ncols,roiIdx);

    x    = prfROIData.(roiName{roiIdx}).eccentricity_smoothed;
    y    = prfROIData.(roiName{roiIdx}).sigma_smoothed;
    w    = prfROIData.(roiName{roiIdx}).varexplained;
    xfit = linspace(opt.mri.eccThresh(1),opt.mri.eccThresh(2),20)';
    yfit = mprf_fit(x,y,w,xfit);
    
    scatter(x, y, [], roi_colors(roiIdx,:),'*'); hold on;
    plot(xfit,yfit,'color','k','LineWidth',4); 

    title(roiName{roiIdx}, 'FontSize', 15)
    ylim([0 opt.mri.eccThresh(2)]); xlim([0 opt.mri.eccThresh(2)]); axis square
    ylabel('PRF size (deg)'); xlabel('PRF eccen (deg');
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH5, fullfile(saveDir,'pRF_size_eccentricity_smoothed_roi_wFits'), '-dpng');


%% Sigma histogram
%-----------------------
fH6 = figure(6); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'Sigma');
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
print(fH6, fullfile(saveDir,'sigma'), '-dpng');

% Smoothed sigma histogram
%-----------------------
fH7 = figure(7); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name','Sigma smoothed');
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
print(fH7, fullfile(saveDir,'sigma_smoothed'), '-dpng');

% Beta histogram
%-----------------------
fH8 = figure(8); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name','Beta')
% maxBeta = prctile(prfROIData.allROIs.beta,[opt.mri.betaPrctileThresh(1) opt.mri.betaPrctileThresh(2)]);

for roiIdx = 1:numRois
    subplot(nrows,ncols,roiIdx);
    dataToPlot = prfROIData.(roiName{roiIdx}).beta;
%     betaMask = dataToPlot<maxBeta(2);
    hist(dataToPlot,nbins)
    title(roiName{roiIdx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    box off;
    ylabel('Frequency'); xlabel('PRF beta (a.u.)')
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH8, fullfile(saveDir,'beta'), '-dpng');

% recomputed beta histogram
%-----------------------
fH9 = figure(9); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'Recomputed beta')
for roiIdx = 1:numRois
    subplot(nrows,ncols,roiIdx);
    dataToPlot = prfROIData.(roiName{roiIdx}).recomp_beta;
    hist(dataToPlot,nbins)
    title(roiName{roiIdx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    box off;
    ylabel('Frequency'); xlabel('PRF Recomp beta (a.u.)')
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end
print(fH9, fullfile(saveDir,'recomp_beta'), '-dpng');


fH10 = figure(10); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920/2   1080], 'Name', 'pRF center distribution');
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
print(fH10, fullfile(saveDir,'prf_center_distribution'), '-dpng');



    
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
    
    vw = loadParameterMap(vw,fullfile(dirPth.fmri.saveDataPth_prfMrv,'nifti', 'varexplained.nii.gz'));
    vemask = vw.map{1}>opt.mri.varExplThresh(1);
    
    for ii = 1:length(prfParams)
        
        fprintf('(%s):  Visualizing %s on mrVista surface \n', mfilename, prfParams{ii})
        
        prfParamNifti =  fullfile(dirPth.fmri.saveDataPth_prfMrv,'nifti', sprintf('%s.nii.gz', prfParams{ii}));
        
        data_type = 'Averages';
        vw = viewSet(vw,'curdt',data_type);  
        vw = viewSet(vw,'displaymode','map');
        vw = loadParameterMap(vw,prfParamNifti);
        vw = viewSet(vw,'cothresh',opt.mri.varExplThresh(1));
        
        
        % Add smoothing/ve mask?       
        vw.map{1} = vw.map{1}.*vemask;
        
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
            otherwise
                cAxisLim = [eps opt.mri.eccThresh(2)];
                vw.ui.mapMode = setColormap(vw.ui.mapMode, 'hsvTbCmap');
                vw = viewSet(vw, 'mapwin', cAxisLim);
                vw = viewSet(vw, 'mapclip', cAxisLim);   
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
                
                figure(fH); clf; cmapbar = vw.ui.mapMode.cmap((vw.ui.mapMode.numGrays+1):end,:);
                imagesc(mrmGet(m, 'screenshot')/255); caxis(cAxisLim); colormap(cmapbar); colorbar; axis image; axis off;
                if opt.saveFig print(fH, fullfile(saveDir,sprintf('%s_%s_%s',idName{thisID}, prfParams{ii},viewList{thisView})), '-dpng'); end
                
            end
        end
    end
end


end

