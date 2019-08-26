function makeAllFigures(subjID)

% Load paths with data files for this subject
dirPth = loadPaths(subjID);

% Go back to root
cd(mprf_rootPath)

% Set options
opt = getOpts;

%% Figures

% Figure 1. Time series (1A) and MEG head plot (1B)
makeFigure1(dirPth,opt);

%Figure 2. Position range line plot and headplots for every position range
opt = getOpts('perturbOrigPRFs', 'position');
makeFigure2(dirPth,opt);

%Figure 3. Size range line plot and headplots for every size range
opt = getOpts('perturbOrigPRFs', 'size');
makeFigure3(dirPth,opt);