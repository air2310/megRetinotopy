function [betas, offsets, varexpl] = regressPredictedResponse(data, prediction, varargin)

% check inputs
if nargin < 3
    regressionType = 'NoOffset';
else
    regressionType = varargin{2};
end
    

% Get design matrix
if strcmp(regressionType, 'NoOffset')
    X = prediction; % do not model offset, because we assume that a blank stimulus should give a 0 response across trials.
elseif strcmp(regressionType, 'WithOffset')
    X = [ones(size(prediction)), prediction];
end

% Regress prediction from phase referenced 10 Hz MEG response
[B,~,~,~,stats] = regress(data,X);

if strcmp(regressionType, 'NoOffset')
    betas   = B(1);
    offsets = [];
elseif strcmp(regressionType, 'WithOffset')
    offsets = B(1);
    betas = B(2);
end
    
varexpl = stats(1);


return
