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


predMEGResponse = predSurfaceResponse * gain';

return
