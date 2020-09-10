function makeSupplementFigureS2(opt)
% Function to create the Supplement Figure S2: MEG head plots with
% variance explained values for individual subjects.

% check inputs
if ~exist('opt', 'var') || isempty(opt)
    opt = getOpts;
end


% Define all subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

% Set up figure for all subjects plotted separately
sz = get(0, 'screensize');
fH1 = figure; set(gcf,'Position', sz);

% Set colormap limits
clims = [0 0.5]; % or clims = [0 max(meanVarExpl)];
interpMethod = 'v4'; % choose 'v4' or 'nearest'

for s = 1:length(subjects)
    
    % Get subject name and directories
    subjectID = subjects{s};
    dirPth = loadPaths(subjectID);
    
    % Find the variance explained file
    varexpl = dir(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp','meanVarExpl.mat'));
    
    % check if the modelPredictions are saved in the folder.
    if isempty(varexpl.folder)
        error('Can''t find mean variance explained mat-file')
    end
    
    % Load the variance explained file
    load(fullfile(varexpl.folder,varexpl.name),'meanVarExpl');
    
    % Plot it!
    subplot(2,5,s)
    megPlotMap(meanVarExpl,clims,fH1, 'parula', ...
        'Mean variance explained', [],[], 'interpmethod', interpMethod);
    c = colorbar;
    c.Location = 'eastoutside';
    c.Box = 'off';
    c.TickDirection = 'out';
    c.TickLength = [0.010 0.010];
    c.FontSize = 12;
end

% Save fig
if opt.saveFig
    % Make folder to save figures
    saveSubDir = ['SuppFigure2'];
    saveDir = fullfile(dirPth.finalFig.savePthAverage, saveSubDir);
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    
    fprintf('\n(%s): Saving Supplemental Figure S2 in %s\n',mfilename, saveDir);
    print(fH1, fullfile(saveDir, sprintf('SuppFigureS2_IndividualSubjects_VarExpOrig%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-depsc');
    print(fH1, fullfile(saveDir, sprintf('SuppFigureS2_IndividualSubjects_VarExpOrig%s_%s_%s', opt.fNamePostFix, sensorsToAverage, summaryMetric)), '-dpng');
end

return
