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

sensorsToAverage = 'top10'; % revert to default: 'top10', choose from: 'allPosterior' or 'top10', 'top5','top15', 'top10Positive', 'top10reliable' 
summaryMetric    = []; % revert to default: 'meanVE', choose from: 
opt              = []; % revert to default made by getOpts()

%% Make Figure 3: SSVEF coherence and reliability
figNum = 3;

makeFigureWrapper(subjects{1}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S1
makeFigureWrapper(subjects{8}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S9 and group average
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, true, summaryMetric, opt) % subject S9 and group average

%% Make Figure 4: Variance Explained by MEG encoding model using original estimated pRFs by fMRI
figNum = 4;

makeFigureWrapper(subjects{1}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S1
makeFigureWrapper(subjects{8}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S9
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, true, summaryMetric, opt)  % group average


%% Make Figure 5: Variance Explained by MEG encoding model using rotated pRF positions
figNum = 5;

makeFigureWrapper(subjects{1}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S1
makeFigureWrapper(subjects{8}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S8
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S9
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, true, summaryMetric, opt)  % group average


%% Make Figure 6: Variance Explained by MEG encoding model using scaled pRF sizes
figNum = 6;

makeFigureWrapper(subjects{1}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S1
makeFigureWrapper(subjects{8}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S8
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, false, summaryMetric, opt) % subject S9
makeFigureWrapper(subjects{9}, figNum, sensorsToAverage, true, summaryMetric, opt)  % group average

%% Make Supplemental Figs

makeSupplementFigureS1
makeSupplementFigureS2
makeSupplementFigureS3(sensorsToAverage);
makeSupplementFigureS4([1,9], true);
makeSupplementFigureS5(sensorsToAverage);
makeSupplementFigureS6([1,9], true);

% Using top 5 and top 15 sensors for position ..
makeSupplementFigureS3('top5');
makeSupplementFigureS3('top15');

% .. and size variations
makeSupplementFigureS5('top5');
makeSupplementFigureS5('top15');

% Using top 10 sensors split half reliable sensors and 
% top 10 sensors from original pRF model fit for position ..
makeSupplementFigureS3('top10reliable');
makeSupplementFigureS3('top10orig');

% .. and size variations
makeSupplementFigureS5('top10reliable');
makeSupplementFigureS5('top10orig'); 

% Sensor-wise group average
makeSupplementFigureS9_sensorwiseAverage_varyPosition(sensorsToAverage, summaryMetric, opt)
makeSupplementFigureS9_sensorwiseAverage_varySize(sensorsToAverage, summaryMetric, opt)

% Variance explained maps using group average 3T retinotopy data from NYU
makeSupplementFigureS10