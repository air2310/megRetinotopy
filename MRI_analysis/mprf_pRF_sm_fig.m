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

for p = 1:length(fn)
    
    allroitmp = [];
    for roiIdx = 1:numRoi
        
        thisFieldName = fn{p};
        data = prf.(fn{p});
        
        if strcmp(thisFieldName,'varexplained')
            data = data(roiLoc{roiIdx});
            prfROIData.(roiName{roiIdx}).(fn{p}) = data(data > opt.mri.varExplThresh(1) & data < opt.mri.varExplThresh(2));
            
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
    prfROIData.allROI.(thisFieldName) = allroitmp;
end

roiName(end+1) = {'allROI'};
numRois = length(roiName);

%% ----------------------
%  1. Visualize distribution of prf parameters for all voxels (whole brain)
%  ----------------------

% Variance explained (all voxels, whole brain)
figure(1); clf; set(gcf, 'Color', 'w', 'Position', [63,894,586,417]);
hist(prf.varexplained,100);
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0.5 0.5];
h.EdgeColor = 'w';
title('Variance explained of pRF model fits (all vertices, whole brain)')
xlabel('Variance explained'); ylabel('Frequency')
set(gca, 'FontSize', 15, 'TickDir', 'out'); box off;

% prf size vs eccentricity (all voxels, whole brain)
figure(2); clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920/3   1080/2]);
smoothing_mask = prf.varexplained > 0;
c_orig = smoothing_mask;
scatter(prf.eccentricity,prf.sigma,[],c_orig);
title('Voxels used for smoothing (ve>0)');
xlabel('Eccentricity  (deg)');
ylabel('Prf size (deg)');
xlim([0 17]); ylim([0 17]); axis square;
set(gca, 'FontSize', 15, 'TickDir', 'out'); box off;
legend('All voxels', 'Included for smoothing (ve>0)'); legend boxoff

%% ------------------------------------------------
%  2. Visualize prf params for different ROIs
%  ------------------------------------------------

nrows = 2;
ncols = ceil(numRois/nrows);
nbins = 50;

% Plot pRF size vs ecc
%-----------------------
figure(3), clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'pRF size vs ecc');
c_orig 	 = [0.5 1 0];
c_smooth = [0 0.5 1];

for roi_idx = 1:numRois
    subplot(nrows,ncols,roi_idx);
    scatter(prfROIData.(roiName{roi_idx}).eccentricity, prfROIData.(roiName{roi_idx}).sigma,[],c_orig,'*');  hold on;
    scatter(prfROIData.(roiName{roi_idx}).eccentricity_smoothed, prfROIData.(roiName{roi_idx}).sigma_smoothed,[],c_smooth,'*');
    
    title(roiName{roi_idx}, 'FontSize', 15)
    legend('original','smoothed','Location','NorthWest')
    legend('Location','NorthWest')
    ylim([0 opt.mri.eccThresh(2)]); xlim([0 opt.mri.eccThresh(2)]); axis square
    ylabel('PRF size (deg)'); xlabel('PRF eccen (deg');
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end


% Sigma histogram
%-----------------------
figure(4), clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'Sigma');
for roi_idx = 1:numRois
    subplot(nrows,ncols,roi_idx);
    hist(prfROIData.(roiName{roi_idx}).sigma,nbins)
    
    title(roiName{roi_idx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    xlim([0 opt.mri.eccThresh(2)]); box off;
    ylabel('Frequency'); xlabel('PRF size (deg)')
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end

% Smoothed sigma histogram
%-----------------------
figure(5), clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name','Sigma smoothed');
for roi_idx = 1:numRois
    subplot(nrows,ncols,roi_idx);
    
    hist(prfROIData.(roiName{roi_idx}).sigma_smoothed,nbins)
    title(roiName{roi_idx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    xlim([0 opt.mri.eccThresh(2)]); box off;
    ylabel('Frequency'); xlabel('PRF smoothed size (deg)')
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end


% Beta histogram
%-----------------------
figure(6), clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name','Beta')
maxBeta = prctile(prfROIData.allROI.beta,[opt.mri.betaPrctileThresh(1) opt.mri.betaPrctileThresh(2)]);

for roi_idx = 1:numRois
    subplot(nrows,ncols,roi_idx);
    dataToPlot = prfROIData.(roiName{roi_idx}).beta;
    betaMask = dataToPlot<maxBeta(2);
    hist(dataToPlot(betaMask),nbins)
    title(roiName{roi_idx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    box off;
    ylabel('Frequency'); xlabel('PRF beta (a.u.)')
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end

% recomputed beta histogram
%-----------------------
figure(7), clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920   1080], 'Name', 'Recomputed beta')
maxBeta = prctile(prfROIData.allROI.recomp_beta,[opt.mri.betaPrctileThresh(1) opt.mri.betaPrctileThresh(2)]);
for roi_idx = 1:numRois
    subplot(nrows,ncols,roi_idx);
    dataToPlot = prfROIData.(roiName{roi_idx}).recomp_beta;
    betaMask = dataToPlot<maxBeta(2);
    hist(dataToPlot(betaMask),nbins)
    title(roiName{roi_idx})
    h = findobj(gca,'Type','patch');
    h.FaceColor = [0 0.5 0.5];
    h.EdgeColor = 'w';
    box off;
    ylabel('Frequency'); xlabel('PRF Recomp beta (a.u.)')
    set(gca, 'FontSize', 14, 'TickDir', 'out')
end

figure(8), clf; set(gcf, 'Color', 'w', 'Position', [10   10   1920/2   1080], 'Name', 'pRF center distribution');

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

for roi_idx = 1:numRois
    
    subplot(nrows,ncols,roi_idx);
    title(roiName{roi_idx})
    polarPlot(0,params); hold on;
    
    x = prfROIData.(roiName{roi_idx}).x;
    y = prfROIData.(roiName{roi_idx}).y;
    x_smoothed = prfROIData.(roiName{roi_idx}).x_smoothed;
    y_smoothed  = prfROIData.(roiName{roi_idx}).y_smoothed;
    
    plot(x,y,'o','MarkerFaceColor',c_orig,'MarkerSize',2); hold all;
    plot(x_smoothed,y_smoothed,'o','MarkerFaceColor',c_smooth,'MarkerSize',2);
    
end

if opt.saveFig
    for f = 1:8
        print(f, fullfile(saveDir,sprintf('fig%d_prf_sm',f)), '-dpng');
    end
end

    
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
        
        
        vw = viewSet(vw,'displaymode','map');
        vw = loadParameterMap(vw,prfParamNifti);
        
        % Add smoothing/ve mask?
        vw.map{1} = vw.map{1}.*smoothing_mask;
        
        switch prfParams{ii}
            case {'polar_angle', 'polar_angle_smoothed'}
                vw = viewSet(vw, 'mapwin', [eps 2*pi]);
                vw = viewSet(vw, 'mapclip', [eps 2*pi]);
                vw.ui.mapMode = setColormap(vw.ui.mapMode, 'hsvCmap');
                
            case {'eccentricity', 'eccentricity_smoothed'}
                % Set colormap and limits
                vw.ui.mapMode = setColormap(vw.ui.mapMode, 'hsvTbCmap');
                vw = viewSet(vw, 'mapwin', [eps opt.mri.eccThresh(2)]);
                vw = viewSet(vw, 'mapclip', [eps opt.mri.eccThresh(2)]);
                
            case {'beta', 'recomp_beta'}
                if opt.mri.betaPrctileThresh
                    maxBetaCmap = prctile(vw.map{1},[opt.mri.betaPrctileThresh(1) opt.mri.betaPrctileThresh(2)]);
                    vw = viewSet(vw, 'mapwin', [eps maxBetaCmap]);
                    vw = viewSet(vw, 'mapclip', [eps maxBetaCmap]);
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
                imagesc(mrmGet(m, 'screenshot')/255); axis image; axis off;
                if opt.saveFig print(fH, fullfile(saveDir,sprintf('%s_%s_%s',idName{thisID}, prfParams{ii},viewList{thisView})), '-dpng'); end
                
            end
        end
    end
end


end

