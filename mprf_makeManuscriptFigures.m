% mprf_makeManuscriptFigures.m
%
% Master script to reproduce manuscript figures

% Define all subjects
subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'}; %

sensorsToAverage = []; % revert to default: 'top10'
summaryMetric    = []; % revert to default: 'meanVE'
opt              = []; % revert to default made by getOpts()

%% Make Figure 3: SSVEF coherence and reliability
figNum = 3;

makeFigureWrapper(subjects{1}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S1
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, true, summaryMetric, opt) % subject S9 and group average

%% Make Figure 4: Variance Explained by MEG encoding model using original estimated pRFs by fMRI
figNum = 4;

makeFigureWrapper(subjects{1}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S1
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S9
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, true, summaryMetric, opt)  % group average


%% Make Figure 5: Variance Explained by MEG encoding model using rotated pRF positions
figNum = 5;

makeFigureWrapper(subjects{1}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S1
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S9
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, true, summaryMetric, opt)  % group average


%% Make Figure 6: Variance Explained by MEG encoding model using scaled pRF sizes
figNum = 6;

makeFigureWrapper(subjects{1}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S1
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S9
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, true, summaryMetric, opt)  % group average

%% Make Supplemental Figs

makeSupplementFigureS1
makeSupplementFigureS2
makeSupplementFigureS3;
makeSupplementFigureS4([1,9], true);
makeSupplementFigureS5;
makeSupplementFigureS6([1,9], true);
