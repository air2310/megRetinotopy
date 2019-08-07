function nIter = checkNumberOfIterations(prf, opt)
% Function to check if the dimensions we are looping over are corresponding
% with the prf parameter estimate we are varying


if strcmp(opt.vary.perturbOrigPRFs, 'position')
    assert(size(prf.x_smoothed_vary,2)==length(opt.vary.position))
    nIter = length(opt.vary.position);
elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
    assert(size(prf.sigma_smoothed_vary,2)==length(opt.vary.size))
    nIter = length(opt.vary.size);
elseif strcmp(opt.vary.perturbOrigPRFs, 'scramble')
    assert(size(prf.sigma_smoothed_scramble,2)==opt.vary.nScrambles)
    nIter = opt.vary.nScrambles;
elseif ~opt.vary.perturbOrigPRFs
    nIter = 1;
end