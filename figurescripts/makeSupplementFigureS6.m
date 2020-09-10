function makeSupplementFigureS6(subjectToPlot, plotAverage)
% Function to make Supplemental Figure S6 of manuscript, showing the effect
% of scaling the original pRF sizes (sigma) on a head plot.

subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

% Define plotting params
clim  = [0 50];
interpmethod = 'v4'; % can also be [] for 'v4' --> smooth interpolation

for s = subjectToPlot
    
    % Get subject ID, options and paths
    subjID  = subjects{s};
    opt     = getOpts('perturbOrigPRFs','size');
    dirPth  = loadPaths(subjID);
    
    % Get range of permutations, plotting rows/columns
    range   = opt.vary.size;
    rows    = 3;
    cols    = round(length(range)/rows)+1;
    
    % Load variance explained
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder, 'pred_resp', 'meanVarExpl'), 'meanVarExpl');
    
    fH1 = figure(1); set(gcf, 'Color', 'w', 'Position', [ 136, 96, 2000,  1138],  'Name', 'Vary pRF size'); hold all;
    
    for ii = 1:length(range)
        % Select data
        meshDataToPlot = meanVarExpl(ii,:).*100;
        
        % Get subplot
        subplot(rows, cols, ii);
        
        megPlotMap(meshDataToPlot,clim,fH1,'parula',...
            sprintf('%1.2fx',range(ii)),[],[],'interpmethod',interpmethod);
        c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
        pos = c.Position; set(c, 'Position', [pos(1)+0.03 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
    end
    
    %% Save figures if requestsed
    if opt.saveFig
        
        [pth, ~] = fileparts(dirPth.model.saveFigPth);
        saveSubDir = ['SupplFigureS6_' opt.subfolder];
        saveDir = fullfile(pth, 'finalfig', saveSubDir);
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end
        
        fprintf('\n(%s): Saving Supplemental Figure S6 in %s\n',mfilename, saveDir);
        
        figure(fH1);
        figurewrite(fullfile(saveDir, sprintf('SupplFigureS6_%s_varySizeMeshes%s', dirPth.subjID, opt.fNamePostFix)),[],[1 300],'.',1);
        %     figurewrite(fullfile(saveDir, sprintf('SupplFigureS6_%s_varySizeMeshes%s', dirPth.subjID, opt.fNamePostFix)),[],0,'.',1);
        
    end
end

if plotAverage
    % Do the same for group average modelfit
    opt = getOpts('perturbOrigPRFs','size');
    plotpRFPerturbationMeshesGroupAverageFit(dirPth, opt)
end
return