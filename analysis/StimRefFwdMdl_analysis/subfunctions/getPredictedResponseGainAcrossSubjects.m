function [predRespSurface, predictedResponseMedianPerROI] = ...
    getPredictedResponseGainAcrossSubjects(subjID, opt)

subjects = {'wlsubj004','wlsubj039', 'wlsubj040', 'wlsubj058','wlsubj068', ...
    'wlsubj070', 'wlsubj081', 'wlsubj106', 'wlsubj109', 'wlsubj111'}; %

sub_dir = strsplit(opt.subfolder,'/');

for s = 1:length(subjects)
   
    dirPth = loadPaths(subjects{s});
      
    predictedResponseDir = fullfile(dirPth.model.saveDataPth, sub_dir{1}, ...
        sub_dir{3}, sub_dir{4});
    
    
    predRespCortex = load(fullfile(predictedResponseDir,'pred_resp', 'predSurfResponseFromPRFs.mat'));
    
    % Load Wang atlas
    roifname      = fullfile(dirPth.fmri.saveDataPth_prfFS, 'pial.wang2015_atlas');
    roiFSWang     = read_curv(roifname);
    areas         = unique(roiFSWang);
    roiFSWang_mask = roiFSWang>0;
    
    predRespCortex_max = max(predRespCortex.predSurfResponse(:,roiFSWang_mask), [],'omitnan');

    for ii = areas(2:end)'
        gainData(s,ii) = median(predRespCortex_max(roiFSWang(roiFSWang_mask)==ii),'omitnan');
    end
   
end

predictedResponseMedianPerROI = mean(gainData,1, 'omitnan');

% Project values onto cortical surface
dirPth = loadPaths(subjID);
predictedResponseDir = fullfile(dirPth.model.saveDataPth, sub_dir{1}, ...
        sub_dir{3}, sub_dir{4});
tmp = load(fullfile(predictedResponseDir,'pred_resp', 'predSurfResponseFromPRFs.mat'));
predRespSurface = NaN(1,size(tmp.predSurfResponse,2));
roifname        = fullfile(dirPth.fmri.saveDataPth_prfFS, 'pial.wang2015_atlas');
roiFSWang       = read_curv(roifname);
    
for ii = areas(2:end)'
    idx = (roiFSWang==ii);
    predRespSurface(idx) = predictedResponseMedianPerROI(ii);
end

end