function makeAllFigures(subjID, whichFigure, sensorsToAverage, plotAverage, summaryMetric, opt)
% Wrapper function for make manuscript figures for MEG Retinotopy project.
% Allows to plot all figures at once, or just a single one, for any
% subject, or the average across subjects
%
% INPUTS:
%   subjID              : subject ID (string)
%   whichFigure         : figures to plot (int or vect), choose from 1,2,3,
%                         or 1:3 (default)
%   sensorsToAverage    : what sensors to average when plotting variations
%                         in size/position. 
%                         Choose from 'allPosterior' (default) or 'top 10'
%   plotAverage         : plot average across subjects (bool), default is
%                         false. When true, individual subject plotting
%                         will be ignored, so any subjectID can be used to
%                         get the opts.
%   summaryMetric       : what summary metric to use when averaging across
%                         subjects (string) Choose from 'meanVE' (default),
%                         'percentChangeVE', 'zscoreVE'.
%
% Example 1: Plot average variance explained for every pRF size variation  
%   makeAllFigures('wlsubj004', 2, 'top10', true, 'meanVE')
% Example 2: Plot average variance explained for every pRF position variation  
%   makeAllFigures('wlsubj004', 1, 'top10', false, 'meanVE')
% Example 3: Plot avergae zcored variance explained for every position pRF variation  
%   makeAllFigures('wlsubj004', 1, 'top10', false, 'zscoreVE')
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019

% Load paths with data files for this subject
dirPth = loadPaths(subjID);

% See what figures to plot
if ~exist('whichFigure','var') || isempty(whichFigure)
    whichFigure = 1:3;
end

if ~exist('sensorsToAverage','var') || isempty(sensorsToAverage)
    sensorsToAverage = 'allPosterior';
end

if ~exist('plotAverage', 'var') || isempty(plotAverage)
    plotAverage = 0;
end

if ~exist('summaryMetric', 'var') || isempty(summaryMetric)
    summaryMetric = 'meanVE';
end

if ~exist('opt', 'var') || isempty(opt)
    
    fullSizeMesh   = true;
    headmodel      = 'OS';
    verbose        = true;
    saveFig        = true;
    addOffsetParam = true;  % if true, use both gain and offset parameters in fit, if false, regress with only 1 free param (gain)
    refitGainParam = true;
    
    % Set options
    opt = getOpts('verbose', verbose, 'saveFig', saveFig, 'headmodel', headmodel, ...
                'fullSizeMesh',fullSizeMesh, 'addOffsetParam', addOffsetParam, 'refitGainParam', refitGainParam);
else
    fullSizeMesh   = opt.fullSizeMesh;
    headmodel      = opt.headmodel;
    verbose        = opt.verbose;
    saveFig        = opt.saveFig;
    addOffsetParam = opt.addOffsetParam;
    refitGainParam = opt.refitGainParam;
    
end

% Go back to root
cd(mprf_rootPath)




%% Figure 1. Time series (1A) and MEG head plot (1B)
if any(intersect(whichFigure,1))
    
    if plotAverage    
        % Average subjects predictions and data separately before fitting
        mprf_CompareGroupAverageDataVsPrediction(dirPth, opt)
        
        % Plot MEG headplots for 10 subjects separately
        makeFigure1Supplement(dirPth, opt) 

    else
    
        % Figure 1. Time series (1A) and MEG head plot (1B)
        makeFigure1(dirPth, opt);
    end
end

%% Figure 2. Position range line plot and headplots for every position range
if any(intersect(whichFigure,2))
    
    opt = getOpts('verbose', verbose, 'saveFig', saveFig, ...
        'perturbOrigPRFs', 'position', 'headmodel', headmodel, ...
        'fullSizeMesh',fullSizeMesh, 'addOffsetParam', addOffsetParam, ...
        'refitGainParam', refitGainParam);
    
    if plotAverage
        makeFigure2AverageSubject(dirPth, opt, sensorsToAverage, summaryMetric);
    else
        makeFigure2(dirPth, opt, sensorsToAverage);
    end
end

%% Figure 3. Size range line plot and headplots for every size range
if any(intersect(whichFigure,3))
    
    opt = getOpts('verbose', verbose, 'saveFig', saveFig, ...
        'perturbOrigPRFs', 'size', 'headmodel', headmodel, ...
        'fullSizeMesh',fullSizeMesh, 'addOffsetParam', addOffsetParam, ...
         'refitGainParam', refitGainParam);
    
    if plotAverage
        makeFigure3AverageSubject(dirPth,opt,sensorsToAverage, summaryMetric);
    else
        makeFigure3(dirPth,opt,sensorsToAverage);
    end
end

