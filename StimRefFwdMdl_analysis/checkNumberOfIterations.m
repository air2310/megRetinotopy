function iter = checkNumberOfIterations(data, opt, type)
% Function to check if the dimensions we are looping over are corresponding
% with the prf parameter estimate we are varying


switch type
    case 'prfSurf' % Prediction on the surface
        
        if strcmp(opt.vary.perturbOrigPRFs, 'position')
            % Check if size data and nr of variations match
            assert(size(data.x_smoothed_vary,2)==length(opt.vary.position))
            
            % Find the original pRFs and place that iteration first, so we
            % can select the vertices after applying the variance threshold
            % to the original data (if we don't do that, and apply the
            % variance threshold to every perturbation of the pRF, the
            % selected vertices will vary between iterations (not good)
            % and sometimes result in no vertices passing the threshold 
            % for example, very small pRFs will have a very large response)
            origIter = (opt.vary.position==0);
            varyIter = 1:length(opt.vary.size);
            varyIter(origIter) = [];          
            iter = [origIter, varyIter]; 
            
        elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
            assert(size(data.sigma_smoothed_vary,2)==length(opt.vary.size))
            origIter = find(opt.vary.size==1);
            varyIter = 1:length(opt.vary.size);
            varyIter(origIter) = [];
            iter = [origIter, varyIter]; 
            
        elseif strcmp(opt.vary.perturbOrigPRFs, 'scramble')
            assert(size(data.sigma_smoothed_scramble,2)==opt.vary.nScrambles)
            iter = 1:opt.vary.nScrambles;
        elseif ~opt.vary.perturbOrigPRFs
            iter = 1;
        end
        
    case {'prfSensor','MEGPhaseRef'} % Prediction for the MEG sensors (after multiplying with the gain matrix)

          % If perturb original pRFs, check dimensions with loaded pRF data
          if strcmp(opt.vary.perturbOrigPRFs, 'position')
              assert(size(data,3)==length(opt.vary.position))
              origIter = (opt.vary.position==0);
              varyIter = 1:length(opt.vary.size);
              varyIter(origIter) = [];          
              iter = [origIter, varyIter]; 
          elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
              assert(size(data,3)==length(opt.vary.size))
              origIter = find(opt.vary.size==1);
              varyIter = 1:length(opt.vary.size);
              varyIter(origIter) = [];
              iter = [origIter, varyIter]; 
          elseif strcmp(opt.vary.perturbOrigPRFs, 'scramble')
              assert(size(predSurfaceResponse,3)==opt.vary.nScrambles)
              iter = 1:opt.vary.nScrambles;
          elseif ~opt.vary.perturbOrigPRFs
              iter = 1;
          end
        
    case 'prfMEGPredvsData'
        phRefAmp10Hz    = data{1};
        predMEGResponse = data{2};
        
        if strcmp(opt.vary.perturbOrigPRFs, 'position')
            assert(size(phRefAmp10Hz,4)==length(opt.vary.position))
            assert(size(predMEGResponse,3)==length(opt.vary.position))
            iter = length(opt.vary.position);
        elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
            assert(size(phRefAmp10Hz,4)==length(opt.vary.size))
            assert(size(predMEGResponse,3)==length(opt.vary.size))
            iter = length(opt.vary.size);
        elseif strcmp(opt.vary.perturbOrigPRFs, 'scramble')
            assert(size(phRefAmp10Hz,4)==opt.vary.nScrambles)
            assert(size(predMEGResponse,3)==opt.vary.nScrambles)
            iter = opt.vary.nScrambles;
        elseif ~opt.vary.perturbOrigPRFs
            iter = 1;
        end
        
end