function [B, varexpl] = regressPredictedResponse(data, prediction)

% Get design matrix (add column of ones to subtract out the mean)
X = [ones(size(prediction)) prediction];
            
% Regress prediction from phase referenced 10 Hz MEG response
[B,~,~,~,stats] = regress(data,X);
            
varexpl = stats(1);

return
