function [] = makeFigure4_SSVEFReliability()
% Function to plot individual subjects for Figure 3C from manuscript,
% plotting 10 Hz SSVEF amplitude split-half reliability, i.e. mean 
% correlation (rho) across 10000 iterations of splitting individual 
% subjects runs randomly into two groups.

subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};


fH2 = figure(2); clf; set(gcf,'Position',[0,300,500,500]); set(fH2, 'Name', 'SSVEF reliability Group Average' , 'NumberTitle', 'off');


saveSubDir = 'Figure4_SSVEFReliability';
saveDir = fullfile(mprf_rootPath,'data','Retinotopy','Quality_check','average','finalfig',saveSubDir);
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

interpMethod = 'v4'; % or if not interpolated, use 'nearest'

% Define opts
opt = getOpts('saveFig', true,'verbose', true, 'fullSizeMesh', true, ...
        'perturbOrigPRFs', false, 'addOffsetParam', false, ...
        'refitGainParam', false);
    
% Get average subject data directory
dirPth = loadPaths(subjects{1});
[dataDir, tmp] = fileparts(dirPth.finalFig.savePthAverage);

% Load all subjects data at once
load(fullfile(dataDir, 'splitHalfAmpReliability1000.mat'), 'splitHalfAmpCorrelation');

%% Plot split half amplitude reliability 
figure(fH1); clf; set(fH1,'Position', [1000, 651, 1500, 687], 'Name', 'SSVEF reliability' , 'NumberTitle', 'off');

for s = 1:length(subjects)  
    subplot(2,5,s);
    ttl = sprintf('S%d', s);
    megPlotMap(splitHalfAmpCorrelation(s,:),[0 .8],fH1, 'hot', ...
        ttl, [],[], 'interpmethod', interpMethod); hold on;
    c = colorbar;c.TickDirection = 'out'; c.Box = 'off';
    pos = c.Position; set(c, 'Position', [pos(1)+0.04 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
%         figurewrite(fullfile(saveDir, sprintf('Figure4_S%d_SSVEFReliabilityCorrelation_All_%s_%s', s, opt.fNamePostFix,interpMethod)),[],0,'.',1);
end

if opt.saveFig
    figure(fH1)
    print(fH1, fullfile(saveDir, sprintf('Figure4_SSVEFReliabilityCorrelation_%s_%s', opt.fNamePostFix,interpMethod)), '-dpng');
    print(fH1, fullfile(saveDir, sprintf('Figure4_SSVEFReliabilityCorrelation_%s_%s', opt.fNamePostFix,interpMethod)),'-pdf');
    fprintf('(%s) Saving figure 3C - All subjects SSVEF reliability in %s\n',mfilename, saveDir);
end

% Plot average split half amplitude reliability
figure(fH2); clf;
megPlotMap(nanmean(splitHalfAmpCorrelation,1),[0 0.8],fH2, 'hot', ...
    'Group Average SSVEF Reliability', [],[], 'interpmethod', interpMethod);
c = colorbar; c.Location='eastoutside';

if opt.saveFig
    print(fH2, fullfile(saveDir, sprintf('Figure4_Group_SSVEFReliabilityCorrelation_%s_%s', opt.fNamePostFix,interpMethod)), '-dpng');
    figurewrite(fullfile(saveDir, sprintf('Figure4_Group_SSVEFReliabilityCorrelation_%s_%s', opt.fNamePostFix,interpMethod)),[],0,'.',1);
    fprintf('(%s): Saving figure 3C - Group Average SSVEF reliability in %s\n',mfilename, saveDir);
end

