function makeSupplementFigureS7(saveFig)
% Function to create the Supplementary Figure S7: 
% Variance explained maps for individual subjects, using a group average
% pRF model from the aggregate NYU3T retinotopy dataset (with 44 subjects).

if nargin<1
    saveFig = 0;
end

% Define all subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

% Set up figure for all subjects plotted separately
% sz = get(0, 'screensize');
% fH1 = figure(1); set(gcf,'Position', sz);

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
    veMEAN(s,:) = meanVarExpl; clear meanVarExpl;
    
    % Now for the original subject specific pRF model
    opt = getOpts('saveFig',saveFig);
     
    % Find the variance explained file
    varexpl = dir(fullfile(dirPth.model.saveDataPth, opt.subfolder,'pred_resp','meanVarExpl.mat'));
    
    % Load the variance explained file
    load(fullfile(varexpl.folder,varexpl.name),'meanVarExpl');
    veSUB(s,:) = meanVarExpl; clear meanVarExpl;
    
    fH1 = figure(fH1); hold on;
    subplot(2,5,s);
    megPlotMap(veMEAN(s,:),clims,fH1, 'parula', ...
        sprintf('S%d', s), [],[], 'interpmethod', interpMethod);
    c = colorbar;
    c.Location = 'eastoutside';
    c.Box = 'off';
    c.Ticks = cticks;
    c.TickDirection = 'out';
    c.TickLength = [0.010 0.010];
    c.FontSize = 12;
    
   saveSubDir = ['SupplFigureS10_' opt.subfolder];
    dirPth = loadPaths('wlsubj004');
    saveDir = fullfile(dirPth.finalFig.savePthAverage, saveSubDir);
    
end


%% Plot histogram of correlation values

veMEAN(isnan(veMEAN))=0;
veSUB(isnan(veSUB))=0;

corrVEmaps = corr(veSUB',veMEAN');
sameSub = diag(corrVEmaps);
idx = logical(1-eye(10));
diffSub = corrVEmaps(idx);


fH3 = figure; h1 = histogram(diffSub, 'Normalization','probability','BinWidth',0.1); hold on;
h2 = histogram(sameSub, 'Normalization','probability','BinWidth',0.1);
h1.FaceColor = [0 0 0];
h2.FaceColor = [1 1 1];
h2.EdgeColor = [0 0 0];
xlim([-0.4, 1]); ylim([0 0.5])
box off; set(gca,'TickDir','out','FontSize',15); axis square
legend({'Across subjects', 'Within subject'})
legend boxoff;
xlabel({'Correlation variance explained maps','using subject-specific vs NYU3T group average pRFs'})
ylabel('Probability')

% Save figures if requested
if saveFig
    saveSubDir = ['SupplFigureS7_' opt.subfolder];
    dirPth = loadPaths('wlsubj004');
    saveDir = fullfile(dirPth.finalFig.savePthAverage, saveSubDir);
    
    fName = sprintf('SupplFigureS7_IndividualSubjects_VarExplOrig_NYU3TAveMaps%s_%s', ...
        opt.fNamePostFix, interpMethod);
    figure(fH1)
    figurewrite(fullfile(saveDir, fName),[],0,'.',1); % EPS
    figurewrite(fullfile(saveDir, fName),[], [1 300],'.',1); % PNG

    fName = sprintf('SupplFigureS7_IndividualSubjects_CorrOrig_vs_NYU3TAveMaps%s', ...
        opt.fNamePostFix);
    figure(fH3)
    figurewrite(fullfile(saveDir, fName),[],0,'.',1); % EPS
    figurewrite(fullfile(saveDir, fName),[], [1 300],'.',1); % PNG
end