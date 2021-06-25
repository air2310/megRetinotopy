% mprf_makeManuscriptFigures.m
%
% Master script to reproduce manuscript figures of
% A population receptive field model of the magnetoencephalography response
% by Kupers, Edadan, Benson, Zuiderbaan, de Jong, Dumoulin and Winawer.
% YEAR. JOURNAL. DOI.
%
% To run all analyses, please see mprf_runAllAnalyses.m
%
% Written by EK @ NYU

%% Define all subjects
subjects = {'wlsubj004','wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'}; %

sensorsToAverage = 'top10'; % [] will revert to default: 'top10' 
                            % choose from: 'allPosterior', 'top5','top10', 
                            % 'top15','top10Positive','top10reliable' 
summaryMetric    = [];      % revert to default: 'meanVE', choose from: 
opt              = [];      % revert to default made by getOpts()

%% Make Figure 3: SSVEF coherence and reliability
figNum = 3;
makeFigureWrapper(subjects{1}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S1
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S9 and group average
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, true, summaryMetric, opt) % subject S9 and group average

%% Make Figure 4: Variance Explained by MEG encoding model using original estimated pRFs by fMRI
figNum = 4;
makeFigureWrapper(subjects{1}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S1
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S9
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, true, summaryMetric, opt)  % group average


%% Make Figure 5: Variance Explained by MEG encoding model using rotated pRF positions
figNum = 5;
sensorSelections = {'top5','top10','top15'};
for senSel = 1:length(sensorSelections)
    makeFigureWrapper(subjects{1}, figNum, sensorSelections{senSel}, false, summaryMetric, opt) % subject S1
    makeFigureWrapper(subjects{9}, figNum, sensorSelections{senSel}, false, summaryMetric, opt) % subject S9
    makeFigureWrapper(subjects{9}, figNum, sensorSelections{senSel}, true, summaryMetric, opt)  % group average
end

%% Make Figure 6: Variance Explained by MEG encoding model using scaled pRF sizes
figNum = 6;
sensorSelections = {'top5','top10','top15'};
for senSel = 1:length(sensorSelections)
    makeFigureWrapper(subjects{1}, figNum, sensorSelections{senSel}, false, summaryMetric, opt)   % subject S1
    makeFigureWrapper(subjects{9}, figNum, sensorSelections{senSel}, false, summaryMetric, opt)   % subject S9
    makeFigureWrapper(subjects{9}, figNum, sensorSelections{senSel}, true, summaryMetric, opt)    % group average modelfit
end

%% Make Figure 7: Fit-then-average group results
makeFigure7A('top10reliable', summaryMetric, opt)
makeFigure7B('top10reliable', summaryMetric, opt)
makeFigure7A('top10', summaryMetric, opt)
makeFigure7B('top10', summaryMetric, opt)

%% Make Supplemental Figs

% S1: Individual subjects - 10 Hz SSVEF coherence and splithalf reliability
makeSupplementFigureS1
% S2: Individual subjects - variance explained by modelfit
makeSupplementFigureS2

% S3: Individual subjects - variance explained by modelfit 
% for pRF position/size variations 
sensorSelections = {'top5','top10','top15'};
for senSel = 1:length(sensorSelections)
    % Using top 5, 10 and 15 sensors for position ..
    makeSupplementFigureS3A(sensorSelections{senSel});
    % .. and size variations
    makeSupplementFigureS3B(sensorSelections{senSel});
end

% S4 & S5: Subject S1, S9 and group variance explained by modelfit
% as spatial topography for pRF position/size variations
makeSupplementFigureS4([1,9], true);
makeSupplementFigureS5([1,9], true);

% S6: Using top 10 sensors split half reliable sensors
% for position ..
makeSupplementFigureS3A('top10reliable');
% .. and size variations
makeSupplementFigureS3B('top10reliable');

% S7: Variance explained maps using a different group average (aggregate 3T
% retinotopy dataset)
makeSupplementFigureS7(false);