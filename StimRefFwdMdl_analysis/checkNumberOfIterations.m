function nIter = checkNumberOfIterations(data, opt, type)
% Function to check if the dimensions we are looping over are corresponding
% with the prf parameter estimate we are varying


switch type
    case 'prfSurf'
        
        if strcmp(opt.vary.perturbOrigPRFs, 'position')
            assert(size(data.x_smoothed_vary,2)==length(opt.vary.position))
            nIter = length(opt.vary.position);
        elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
            assert(size(data.sigma_smoothed_vary,2)==length(opt.vary.size))
            nIter = length(opt.vary.size);
        elseif strcmp(opt.vary.perturbOrigPRFs, 'scramble')
            assert(size(data.sigma_smoothed_scramble,2)==opt.vary.nScrambles)
            nIter = opt.vary.nScrambles;
        elseif ~opt.vary.perturbOrigPRFs
            nIter = 1;
        end
        
    case 'MEGPhaseRef'
        if strcmp(opt.vary.perturbOrigPRFs, 'position')
            assert(size(data,3)==length(opt.vary.position))
            nIter = length(opt.vary.position);
        elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
            assert(size(data,3)==length(opt.vary.size))
            nIter = length(opt.vary.size);
        elseif strcmp(opt.vary.perturbOrigPRFs, 'scramble')
            assert(size(data,3)==opt.vary.nScrambles)
            nIter = opt.vary.nScrambles;
        elseif ~opt.vary.perturbOrigPRFs
            nIter = 1;
        end
        
    case 'prfMEG'
        phRefAmp10Hz    = data{1};
        predMEGResponse = data{2};
        
        if strcmp(opt.vary.perturbOrigPRFs, 'position')
            assert(size(phRefAmp10Hz,4)==length(opt.vary.position))
            assert(size(predMEGResponse,3)==length(opt.vary.position))
            nIter = length(opt.vary.position);
        elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
            assert(size(phRefAmp10Hz,4)==length(opt.vary.size))
            assert(size(predMEGResponse,3)==length(opt.vary.size))
            nIter = length(opt.vary.size);
        elseif strcmp(opt.vary.perturbOrigPRFs, 'scramble')
            assert(size(phRefAmp10Hz,4)==opt.vary.nScrambles)
            assert(size(predMEGResponse,3)==opt.vary.nScrambles)
            nIter = opt.vary.nScrambles;
        elseif ~opt.vary.perturbOrigPRFs
            nIter = 1;
        end
        
end