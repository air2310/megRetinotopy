function makeSupplementFigureS4(subjectToPlot, plotAverage)
% Function to visualize variance explained meshes when varying pRF position
% for the request subject and group average fit
%
%
% Written by Eline Kupers, NYU 2020

subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

for s = subjectToPlot
    
    % Get subject ID, options and paths
    subjID  = subjects{s};
    opt     = getOpts('perturbOrigPRFs','size');
    dirPth  = loadPaths(subjID);
    
    % Define the range of rotations
    range   = opt.vary.position;
    
    % Define figure/subfigure locations
    rows = 2;
    cols = round(length(range)/rows);
    fH1 = figure(1); clf; set(fH1, 'Color', 'w', 'Position', [326,584,1234,754], 'Name', 'Vary pRF position'); hold all;
    
    % plotting params
    clim  = [0 50];
    interpmethod = 'v4'; % can be [] for 'v4' --> smooth interpolation or 'nearest' to avoid interpolation
    
    
    for ii = 1:length(range)
        
        % Select data
        meshDataToPlot = meanVarExpl(ii,:).*100;
        
        % Get subplot
        hold all;
        subplot(rows, cols, ii);
        
        megPlotMap(meshDataToPlot,clim,fH1,'parula',...
            range(ii),[],[],'interpmethod',interpmethod);
        c = colorbar; c.TickDirection = 'out'; c.Box = 'off';
        pos = c.Position; set(c, 'Position', [pos(1)+0.04 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
        
    end
    
    % Save fig if requested
    if opt.saveFig
        
        [pth, ~] = fileparts(dirPth.model.saveFigPth);
        saveSubDir = ['SupplFigureS4_' opt.subfolder];
        saveDir = fullfile(pth, 'finalfig', saveSubDir);
        if ~exist(saveDir, 'dir')
            mkdir(saveDir);
        end
        
        fprintf('\n(%s): Saving Supplemental Figure S4 in %s\n',mfilename, saveDir);
        
        figurewrite(fullfile(saveDir, sprintf('SupplFigureS4_%s_varyPositionMeshes%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)),[],[1 300],'.',1);
        figurewrite(fullfile(saveDir, sprintf('SupplFigureS4_%s_varyPositionMeshes%s_%s', dirPth.subjID, opt.fNamePostFix, sensorsToAverage)),[],0,'.',1);
    end
end

if plotAverage
    % Do the same for group average modelfit
    opt = getOpts('perturbOrigPRFs','position');
    plotpRFPerturbationMeshesGroupAverageFit(dirPth, opt)
    
end
return