function [] = makeFigure4_SSVEFReliability()



subjects = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'};


fH1 = figure(1); clf; set(gcf,'Position',[0,300,500,500]); set(fH1, 'Name', 'SSVEF reliability' , 'NumberTitle', 'off');
fH2 = figure(2); clf; set(gcf,'Position',[0,300,500,500]); set(fH1, 'Name', 'SSVEF reliability Group Average' , 'NumberTitle', 'off');


saveSubDir = 'Figure4_SSVEFReliability';
saveDir = fullfile(mprf_rootPath,'data','Retinotopy','Quality_check','average','finalfig',saveSubDir);

interpMethod = 'v4'; % or if not interpolated, use 'nearest'

if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

figure(fH1); clf; set(gcf,'Position',[1000, 651, 1500, 687]);


for s = 1:length(subjects)
    
    % Get subject ID, options and paths
    subjID = subjects{s};
    opt    = getOpts('verbose', true, 'saveFig',true, 'headmodel', 'OS','fullSizeMesh',true);
    dirPth = loadPaths(subjID);
    
    % Load data files for this subject
    load(fullfile(dirPth.model.saveDataPth, opt.subfolder, 'refphase', 'splitHalfAmpReliability.mat'), 'splitHalfAmpReliability');
    allSubjectData(s,:) = splitHalfAmpReliability;
    
    
    % Plot split half amplitude reliability
    subplot(2,5,s);
    ttl = sprintf('S%d', s);
    megPlotMap(splitHalfAmpReliability,[0 .8],fH1, 'hot', ...
        ttl, [],[], 'interpmethod', interpMethod); hold on;
    c = colorbar;c.TickDirection = 'out'; c.Box = 'off';
    pos = c.Position; set(c, 'Position', [pos(1)+0.04 pos(2)+0.03, pos(3)/1.5, pos(4)/1.5])
end    
    if opt.saveFig
        figure(fH1)
        print(fH1, fullfile(saveDir, sprintf('Figure4_S%d_SSVEFReliability_All_%s',  opt.fNamePostFix,interpMethod)), '-dpng');
        figurewrite(fullfile(saveDir, sprintf('Figure4_S%d_SSVEFReliability_All_%s',  opt.fNamePostFix,interpMethod)),[],0,'.',1);
        fprintf('\n saving figure XX: Individual Subject SSVEF reliability in %s',saveDir);
    end
    



% Plot split half amplitude reliability
figure(fH2); clf;
megPlotMap(nanmean(allSubjectData,1),[0 0.5],fH2, 'hot', ...
    'Group Average SSVEF Reliability', [],[], 'interpmethod', interpMethod);
c = colorbar; c.Location='eastoutside';

if opt.saveFig
    print(fH2, fullfile(saveDir, sprintf('Figure4_Group_SSVEFReliability_%s_%s', opt.fNamePostFix,interpMethod)), '-dpng');
    figurewrite(fullfile(saveDir, sprintf('Figure4_Group_SSVEFReliability_%s_%s', opt.fNamePostFix,interpMethod)),[],0,'.',1);
    fprintf('\n saving figure XX: Group Average SSVEF reliability in %s',saveDir);
end

