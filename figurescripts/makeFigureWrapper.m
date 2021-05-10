function makeFigureWrapper(subjID, whichFigure, sensorsToAverage, plotAverage, summaryMetric, opt)
% Wrapper function for make manuscript figures for MEG Retinotopy project.
% Allows to plot all figures at once, or just a single one, for any
% subject, or the average across subjects
%
% INPUTS:
%   subjID              : subject ID (string)
%   whichFigure         : figures to plot (int or vect), choose from 1,2,3,
%                         or 1:3 (default)
%   sensorsToAverage    : what sensors to average when plotting variations
%                         in size/position. Choose from 'allPosterior',
%                         'allPosteriorPositive', 'top10' (default, union
%                         of top 10 sensors across all iterations), 'top5',
%                         'top15', (same computation as top10, but less or
%                         more sensors),'top10Positive' (only positive top
%                         10 sensors at each iteration), or 'top10reliable'
%                         use sensors with highest split half correlation.
%   plotAverage         : plot average across subjects (bool), default is
%                         false. When true, individual subject plotting
%                         will be ignored, so any subjectID can be used to
%                         get the opts.
%   summaryMetric       : what summary metric to use when averaging across
%                         subjects (string) Choose from 'meanVE' (default),
%                         'percentChangeVE', 'zscoreVE'.
%
% Example 1: Plot average variance explained for every pRF size variation  
%   makeFigureWrapper('wlsubj004', 4, 'top10', true, 'meanVE', [])
% Example 2: Plot average variance explained for every pRF position variation  
%   makeFigureWrapper('wlsubj004', 5, 'top10', false, 'meanVE', [])
% Example 3: Plot avergae zcored variance explained for every position pRF variation  
%   makeFigureWrapper('wlsubj004', 6, 'top10', false, 'zscoreVE', [])
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019

% Load paths with data files for this subject
dirPth = loadPaths(subjID);

% See what figures to plot
if ~exist('whichFigure','var') || isempty(whichFigure)
    whichFigure = 3:6;
end

if ~exist('sensorsToAverage','var') || isempty(sensorsToAverage)
    sensorsToAverage = 'top10';
end

if ~exist('plotAverage', 'var') || isempty(plotAverage)
    plotAverage = false;
end

if ~exist('summaryMetric', 'var') || isempty(summaryMetric)
    summaryMetric = 'meanVE';
end

if ~exist('opt', 'var') || isempty(opt)
    
    fullSizeMesh   = true;
    headmodel      = 'OS';
    verbose        = true;
    saveFig        = true;
    addOffsetParam = false;  % if true, use both gain and offset parameters in fit, if false, regress with only 1 free param (gain)
    refitGainParam = false;
    
    % Set options
    opt = getOpts('verbose', verbose, 'saveFig', saveFig, 'headmodel', headmodel, ...
                    'fullSizeMesh',fullSizeMesh, 'addOffsetParam', addOffsetParam, ...
                    'refitGainParam', refitGainParam);
else
    fullSizeMesh   = opt.fullSizeMesh;
    headmodel      = opt.headmodel;
    verbose        = opt.verbose;
    saveFig        = opt.saveFig;
    addOffsetParam = opt.addOffsetParam;
    refitGainParam = opt.refitGainParam;
    useHCPAveMaps  = opt.mri.useHCPAveMaps;
end

% Go back to root
cd(mprf_rootPath)

%% Figure 3. Steady state coherence / reliability
if any(intersect(whichFigure,3))
    
    subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};

    subjectToPlot = find(strcmp(subjects,subjID));
    
    plotSupplementalFig = false;
    
    % Figure 3A. Steady state coherence
    makeFigure3AB_SSVEFCoherence(subjectToPlot,plotAverage,plotSupplementalFig)
    
    % Figure 3B. Steady state reliability
    makeFigure3C_SSVEFReliability(subjectToPlot,plotAverage,plotSupplementalFig)

end



%% Figure 4. Time series (4A) and MEG head plot (4B)
if any(intersect(whichFigure,4))
    
    if plotAverage    
        % Average subjects predictions and data separately before fitting
        plotGroupAverageDataVsPrediction(dirPth, opt, whichFigure, sensorsToAverage)
    else  
        % Figure 4. Time series (4A) and MEG head plot (4B)
        makeFigure4(dirPth, opt);
    end
end

%% Figure 5. Position range line plot and headplots for every position range
if any(intersect(whichFigure,5))
    
    opt = getOpts('verbose', verbose, 'saveFig', saveFig, ...
        'perturbOrigPRFs', 'position', 'headmodel', headmodel, ...
        'fullSizeMesh',fullSizeMesh, 'addOffsetParam', addOffsetParam, ...
        'refitGainParam', refitGainParam, 'useHCPAveMaps', useHCPAveMaps);
    
    if plotAverage
        % Group Average modelfit variance explained (averaging subjects 
        % predictions and data separately before fitting)
        plotGroupAverageDataVsPrediction(dirPth, opt, whichFigure, sensorsToAverage)
    else
        makeFigure5(dirPth, opt, sensorsToAverage);
    end
end

%% Figure 6. Size range line plot and headplots for every size range
if any(intersect(whichFigure,6))
    
    opt = getOpts('verbose', verbose, 'saveFig', saveFig, ...
        'perturbOrigPRFs', 'size', 'headmodel', headmodel, ...
        'fullSizeMesh',fullSizeMesh, 'addOffsetParam', addOffsetParam, ...
         'refitGainParam', refitGainParam, 'useHCPAveMaps', useHCPAveMaps);
    
    if plotAverage       
        % Group Average modelfit (averaging subjects predictions and data
        % separately before fitting)
        plotGroupAverageDataVsPrediction(dirPth, opt, whichFigure, sensorsToAverage)
    else
        makeFigure6(dirPth,opt,sensorsToAverage);
    end
end

