function [prfParams, nIter] = getpRFParamNames(opt)
% Function to get the appropiate prf parameter file names to load..

if strcmp(opt.vary.perturbOrigPRFs, 'position')
    prfParams = {'varexplained', 'mask', 'recomp_beta', 'x_smoothed_vary.mgz', 'y_smoothed_vary.mgz', 'sigma_smoothed'};
elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
    prfParams = {'varexplained', 'mask', 'recomp_beta', 'x_smoothed', 'y_smoothed', 'sigma_smoothed_vary.mgz'};
elseif strcmp(opt.vary.perturbOrigPRFs, 'scramble')
    prfParams = {'varexplained', 'mask', 'recomp_beta_scramble.mgz', 'x_smoothed_scramble.mgz', 'y_smoothed_scramble.mgz', 'sigma_smoothed_scramble.mgz'};
elseif (~opt.mri.useBensonMaps && opt.mri.useSmoothedData)
    prfParams = {'varexplained', 'mask', 'recomp_beta', 'x_smoothed', 'y_smoothed', 'sigma_smoothed'};
elseif (~opt.mri.useBensonMaps && opt.mri.useSmoothedData && opt.roi.onlyV123WangAtlas)
    prfParams = {'varexplained', 'V123mask', 'recomp_beta', 'x_smoothed', 'y_smoothed', 'sigma_smoothed'};
elseif opt.mri.useBensonMaps
    prfParams = {'mask', 'beta', 'x', 'y', 'sigma'};
elseif opt.roi.onlyV123WangAtlas
    prfParams = {'varexplained', 'V123mask', 'x', 'y', 'sigma', 'beta'};
else
    prfParams = {'varexplained', 'mask', 'x', 'y', 'sigma', 'beta'};
end