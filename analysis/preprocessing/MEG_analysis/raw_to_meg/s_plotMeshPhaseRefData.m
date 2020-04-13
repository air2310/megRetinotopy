
subjID = 'wlsubj081';
dirPth = loadPaths(subjID);
type = 'incoherent';
fname = sprintf('/Volumes/server/Projects/MEG/Retinotopy/Quality_check/%s/modelfit/original/%s/pred_resp/phaseReferencedMEGData.mat', subjID, type);
load(fname)

rows = 2;
cols = ceil(size(phRefAmp10Hz,2)/rows);

fH = figure(1); clf; set(fH, 'Position', [17          38        2371        1300]); 
for ii = 1:size(phRefAmp10Hz,2)
    
    subplot(rows, cols, ii)
    
    dataToPlot = squeeze(nanvar(phRefAmp10Hz(:,ii,:)));
    megPlotMap(dataToPlot,[], [],[],[],[],[],'interpmethod', 'nearest');
    title(sprintf('Split %d', ii))
end

print(gcf, '-dpng', fullfile(dirPth.model.saveFigPth,'original', type, sprintf('%s_PhaseReferencedAmps_var', subjID)))
