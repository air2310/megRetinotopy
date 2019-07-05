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

if opt.verbose
    figure, set(gcf, 'Position', [652   784   908   554], 'Color', 'w');
    plot(predMEGResponse); hold on;
    plot(zeros(1,size(predMEGResponse,1)),'k')
    title('Predicted MEG response from MRI prfs')
    xlabel('time points (epochs)');
    ylabel('MEG response (T??)');
    ylim([-1,1]*max(abs(predMEGResponse(:))))
    set(gca, 'FontSize', 14, 'TickDir', 'out'); box off;
end


return
