function predResponseAllVertices = mprf_MEGPredictionFromSurface(prf, stim)
% Function to make predictions for every vertex on given surface
%
%   predResponseAllVertices = mprf_MEGPredictionFromSurface(prf, stim)
%
% INPUTS:
%   prf               : struct with prf data, separate for every parameter
%   stim              : meg stimulus struct (should contain fieldname
%                      'im', with windowed stimulus images for every epoch)
%
% OUTPUTS:
%   predResponseAllVertices : predicted response time series (vertex x timepoints)
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019



% Find loaded prf parameters
fn = fieldnames(prf);
% x0    = prf.(fn{cellfind(regexp(fn, '\<x'))});
% y0    = prf.(fn{cellfind(regexp(fn, '\<y'))});
% sigma = prf.(fn{cellfind(regexp(fn, '\<sigma'))});
% beta  = prf.(fn{cellfind(regexp(fn, 'beta'))});

x0    = prf.(fn{find(~cellfun(@isempty,regexp(fn, '\<x')))});
y0    = prf.(fn{find(~cellfun(@isempty,regexp(fn, '\<y')))});
sigma = prf.(fn{find(~cellfun(@isempty,regexp(fn, '\<sigma')))});
beta  = prf.(fn{find(~cellfun(@isempty,regexp(fn, 'beta')))});

% Preallocate space
predResponseAllVertices = NaN(size(stim.im,2), size(prf.varexplained,1));

% Build RFs from prf surface parameters that fall within the stimulus aperture (stim locations x vertices)
RF = rfGaussian2d(stim.X, stim.Y, sigma, sigma, false, x0, y0);

% Rescale beta's by prf size (sigma) -- since larger sizes will produce
% larger responses
% beta = beta./sqrt(2*pi*sigma.^2);
    
% Get predicted response from RFs given stimulus (epochs x vertices)
predResponse = bsxfun(@times, stim.im' * RF, beta');

% Only output the predicted response for vertices within masks
predResponseAllVertices(:, prf.vemask & prf.roimask & prf.eccenmask) = predResponse;


return
