function iter = checkNumberOfIterations(data, opt, type)
% Function to check if the data dimensions we are looping over. In case we 
% perturb the original prf parameters estimates (i.e. varying the location,
% varying the size or scrambling the vertex locations), we need to check if
% the last data dimension matches the number of pertubations. 
% 
% When computing the predicted time series on the cortical surface, we will
% put the original estimated prf parameters first in the iteration list. We
% do this, because we only want to apply the variance threshold once, on
% the original prf estimates. If we apply this variance threshold to all
% iterations of perturbed prfs, the selected vertices will vary between 
% iterations (not desirable) and/or there is a chance that none of the
% vertex' predictions will survive. For example, very small pRFs will have 
% a very large predicted response and therefore seen as outliers and set to 
% NaNs. This would defeat the purpose of the pertubation analysis.
%   
%   iter = checkNumberOfIterations(data, opt, type)
%
% INPUTS:
%   data    : data to check its dimensions, data type can vary:
%              - predicted surface time series (prfSurf)
%                  (nr vertices x nr iterations)
%              - predicted MEG sensor time series (prfSensor)
%                  (nr epochs x nr sensors x nr iterations)
%              - observed MEG sensor amplitudes (MEGPhaseRef)
%                  (nr epochs x nr sensors x nr iterations)
%              - both predicted and observed MEG sensor amplitudes (prfMEGPredvsData)
%                  (nr epochs x nr sensors x nr iterations)
%   opt     : struct with boolean flag options
%   type    : string to define type of data to check, can be:
%              'prfSurf', 'prfSensor', 'MEGPhaseRef', 'prfMEGPredvsData'
%              see data input variable for array dimensions
%
%
% OUTPUT:
%   iter    : vector of for-loop iteration order prf pertubations (1 x n)
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019

switch type
    case 'prfSurf' % PRF time course predictions on the cortical surface
        
        if strcmp(opt.vary.perturbOrigPRFs, 'position')
            % Check if size data and nr of variations match
            if opt.mri.useSmoothedData
                assert(size(data.x_smoothed_vary,2)==length(opt.vary.position))
            else
                assert(size(data.x_vary,2)==length(opt.vary.position))
            end
            % Find the original pRFs and place that iteration first, so we
            % can select the vertices after applying the variance threshold
            % to the original data.
            origIter = find(opt.vary.position==0);
            varyIter = 1:length(opt.vary.position);
            varyIter(origIter) = [];          
            iter     = [origIter, varyIter]; 
            
        elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
            % We do the same in case of size variations (polar angle rotations)
            % i.e. check size of data, get the iteration nr of original prf
            % params, place that one first in the iteration list, so we can
            % apply our variance threshold.
            if opt.mri.useSmoothedData
                assert(size(data.sigma_smoothed_vary,2)==length(opt.vary.size))
            else
                assert(size(data.sigma_vary,2)==length(opt.vary.size))
            end
            origIter = find(opt.vary.size==1);
            varyIter = 1:length(opt.vary.size);
            varyIter(origIter) = [];
            iter     = [origIter, varyIter]; 
            
        elseif strcmp(opt.vary.perturbOrigPRFs, 'scramble')
            % Again, check size of data, we add 1 iteration to account for  
            % original estimated prf data placed as the first set of
            % vertices
            if opt.mri.useSmoothedData
                assert(size(data.sigma_smoothed_scramble,2)==opt.vary.nScrambles+1)
            else
                assert(size(data.sigma_scramble,2)==opt.vary.nScrambles+1)
            end
            iter = 1:(opt.vary.nScrambles+1); 
        
        elseif ~opt.vary.perturbOrigPRFs
            iter = 1;
        end
        
    case {'prfSensor','MEGPhaseRef'} % Prediction for the MEG sensors (after multiplying with the gain matrix)

          % We do not need to change the order of the for loop iterations,
          % because we only apply the variance threshold to catch outliers
          % at the surface prediction stage.
          if strcmp(opt.vary.perturbOrigPRFs, 'position')
              assert(size(data,3)==length(opt.vary.position))
              iter = 1:length(opt.vary.position);         
          elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
              assert(size(data,3)==length(opt.vary.size))
              iter = 1:length(opt.vary.size);
          elseif strcmp(opt.vary.perturbOrigPRFs, 'scramble')
              assert(size(predSurfaceResponse,3)==opt.vary.nScrambles+1)
              iter = 1:(opt.vary.nScrambles+1);
          elseif ~opt.vary.perturbOrigPRFs
              iter = 1;
          end
        
    case 'prfMEGPredvsData' % When comparing predicted PRF time series (in sensor space) to observed MEG data (in sensor space)
        
        phRefAmp10Hz    = data{1};
        predMEGResponse = data{2};
        
        % We do not need to change the order of the for loop iterations,
        % because we only apply the variance threshold to catch outliers
        % at the surface prediction stage.
        if strcmp(opt.vary.perturbOrigPRFs, 'position')
            assert(size(phRefAmp10Hz,4)==length(opt.vary.position))
            assert(size(predMEGResponse,3)==length(opt.vary.position))
            iter = 1:length(opt.vary.position);
        elseif strcmp(opt.vary.perturbOrigPRFs, 'size')
            assert(size(phRefAmp10Hz,4)==length(opt.vary.size))
            assert(size(predMEGResponse,3)==length(opt.vary.size))
            iter = 1:length(opt.vary.size);
        elseif strcmp(opt.vary.perturbOrigPRFs, 'scramble')
            assert(size(phRefAmp10Hz,4)==opt.vary.nScrambles+1)
            assert(size(predMEGResponse,3)==opt.vary.nScrambles+1)
            iter = 1:(opt.vary.nScrambles+1);
        elseif ~opt.vary.perturbOrigPRFs
            iter = 1;
        end
        
end