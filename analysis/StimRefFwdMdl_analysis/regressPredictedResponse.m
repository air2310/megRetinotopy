function [betas, offsets, varexpl] = regressPredictedResponse(data, prediction, varargin)

% check inputs
addOffsetParam = varargin{2};

% Get design matrix
if addOffsetParam
    X = [ones(size(prediction)), prediction];
else
    X = prediction; % do not model offset, because we assume that a blank stimulus should give a 0 response across trials.
end

% Regress prediction from phase referenced 10 Hz MEG response
[B,~,~,~,stats] = regress(data,X);

if addOffsetParam
    offsets = B(1);
    betas   = B(2);
else
    offsets = [];
    betas   = B(1);
end

varexpl = stats(1);


return
