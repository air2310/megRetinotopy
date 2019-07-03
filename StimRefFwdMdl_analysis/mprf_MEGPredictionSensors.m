function predMEGResponse = mprf_MEGPredictionSensors(predSurfaceResponse, gain)
% Function to compute predicted response at MEG sensors from predicted
% surface response (using MEG stimulus), by multiplying responses with gain
% matrix
% 
%   predMEGResponse = mprf_MEGPredictionSensors(predResponse, gain, subjID) 
%
% INPUTS:
%   predSurfaceResponse :  predicted responses from surface (epochs x vertices)
%   gain                :  gain matrix, weighted sum of vertices contributing 
%                           to each individual MEG sensor (sensors x vertices)
%
% OUTPUTS:
%   predicted response time series for every MEG sensor (epochs x sensors)

% Get nans in matrix
nanMask = isnan(predSurfaceResponse);

% Create output matrix
% predMEGResponse = NaN(size(nanMask,1), size(gain,1));

% Replace NaNs with zero's output matrix
predSurfaceResponse(nanMask)=0;

predMEGResponse = predSurfaceResponse * gain';


return
