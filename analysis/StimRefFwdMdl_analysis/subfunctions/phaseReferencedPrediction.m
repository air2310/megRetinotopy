function [err, pred_out] = phaseReferencedPrediction(ref_phase, pred, A, P, hasOffset)

% Get phase referenced data
phRef10Hz = rescaleAmpsWithRefPhase(A, P, ref_phase);

% Get phase referenced data
[B, offset] = regressPredictedResponse(phRef10Hz, pred, 'addOffsetParam', hasOffset);

% Get new model prediction, scaled with gain and offset
pred_out = pred * B + offset;

% Get error (SSres / SStot)
if isempty(A)
    err = [];
else
    err = sum((pred_out - phRef10Hz).^2) / sum((phRef10Hz - mean(phRef10Hz)).^2);
end

end
