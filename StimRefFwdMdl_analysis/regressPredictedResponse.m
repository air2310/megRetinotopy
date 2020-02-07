function [B, varexpl] = regressPredictedResponse(data, prediction)

% Get design matrix (add column of ones to subtract out the mean)
% X = [ones(size(prediction)) prediction];
X = prediction; % do not model offset, because we assume that a blank stimulus should give a 0 response across trials.


% Regress prediction from phase referenced 10 Hz MEG response
[B,~,~,~,stats] = regress(data,X);
            
varexpl = stats(1);

return
