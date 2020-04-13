function predMEGResponse = mprf_MEGPredictionSensors(predSurfaceResponse, gain)
% Function to compute predicted response at MEG sensors from predicted
% surface response (using MEG stimulus), by multiplying responses with gain
% matrix
% 
%   predMEGResponse = mprf_MEGPredictionSensors(predSurfaceResponse, gain) 
%
% INPUTS:
%   predSurfaceResponse : predicted responses from surface
%                           (epochs x vertices)
%   gain                : gain matrix, weighted sum of vertices contributing 
%                           to each individual MEG sensor
%                           (sensors x vertices)
%
% OUTPUT:
%   predMEGResponse     : predicted MEG sensor responses
%                           (epochs x sensors)
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019

% Locate nans in matrix
nanMask = isnan(predSurfaceResponse);

% Replace NaNs with zero's output matrix
predSurfaceResponse(nanMask)=0;

% Get predict MEG response by multiplying response with gain matrix 
predMEGResponse = predSurfaceResponse * gain';

return
