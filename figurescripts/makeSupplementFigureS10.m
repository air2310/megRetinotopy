function makeSupplementFigureS10(saveFig)
% Function to create the Supplement Figure S10: 
% Variance explained maps for individual subjects, using a group average
% pRF model from the NYU3T retinotopy dataset (with 44 subjects).


% Define all subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

% Set up figure for all subjects plotted separately
sz = get(0, 'screensize');
fH1 = figure(1); set(gcf,'Position', sz);
fH2 = figure(2); set(gcf,'Position', sz);

% Set colormap limits
clims  = [0 0.5]; % or clims = [0 max(meanVarExpl)];
cticks = [0 0.25 0.5];
interpMethod = 'v4'; % choose 'v4' or 'nearest'
interplim    = 'convex';

for s = 1:length(subjects)
    
    % Get subject name and directories
    subjectID = subjects{s};
    dirPth = loadPaths(subjectID);
    
    % Get variance explained by NYU 3T group average pRF model
    opt = getOpts('useNYU3TAveMaps', true, 'saveFig',saveFig);
    
    % Find the variance explained file
    varexpl = dir(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp','meanVarExpl.mat'));
    
    % Load the variance explained file
    load(fullfile(varexpl.folder,varexpl.name),'meanVarExpl');
    ve_nyu3t(s,:) = meanVarExpl; clear meanVarExpl;
    
    % Now for the original subject specific pRF model
    opt = getOpts('saveFig',saveFig);
     
    % Find the variance explained file
    varexpl = dir(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp','meanVarExpl.mat'));
    
    % Load the variance explained file
    load(fullfile(varexpl.folder,varexpl.name),'meanVarExpl');
    ve_orig(s,:) = meanVarExpl; clear meanVarExpl;
    
    figure(fH1); hold on;
    subplot(2,5,s)
    megPlotMap(ve_nyu3t(s,:),clims,fH1, 'parula', ...
        sprintf('S%d', s), [],[], 'interpmethod', interpMethod);
    c = colorbar;
    c.Location = 'eastoutside';
    c.Box = 'off';
    c.Ticks = cticks;
    c.TickDirection = 'out';
    c.TickLength = [0.010 0.010];
    c.FontSize = 12;
end

%% Now plot unique sensors with variance explained > 10%
mask1  = ve_orig>0.1;
mask2  = ve_nyu3t>0.1;
diffMask = mask1-mask2;

clims = [-1,1];
cticks = [-1,0,1];

for s = 1:length(subjects)
    figure(fH2); subplot(2,5,s)
    megPlotMap(diffMask(s,:),clims,fH2, 'bipolar', ...
        sprintf('S%d', s), [],[], 'interpmethod', 'nearest', 'interplim',interplim);
    c = colorbar;
    c.Location = 'eastoutside';
    c.Box = 'off';
    c.Ticks = cticks;
    c.TickDirection = 'out';
    c.TickLength = [0.010 0.010];
    c.FontSize = 12;
end


% Save figures if requested
if saveFig
    saveSubDir = ['SupplFigureS10_' opt.subfolder];
    dirPth = loadPaths('wlsubj004');
    saveDir = fullfile(dirPth.finalFig.savePthAverage, saveSubDir);
    
    fName = sprintf('SupplFigureS10_IndividualSubjects_VarExplOrig_NYU3TAveMaps%s_%s', ...
        opt.fNamePostFix, interpMethod);
    figure(fH1)
    figurewrite(fullfile(saveDir, fName),[],0,'.',1); % EPS
    figurewrite(fullfile(saveDir, fName),[], [1 300],'.',1); % PNG
    
    fName = sprintf('SupplFigureS10_IndividualSubjects_SetDiffOrig_NYU3TAveMaps%s_%s_mask10prct', ...
        opt.fNamePostFix, interpMethod);
    figure(fH2)
    figurewrite(fullfile(saveDir, fName),[],0,'.',1); % EPS
    figurewrite(fullfile(saveDir, fName),[], [1 300],'.',1); % PNG
end