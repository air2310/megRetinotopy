function pred = mprf__load_model_predictions

pred = [];

global mprfSESSION
if isempty(mprfSESSION)
    load('mprfSESSION.mat')
end

try
    pred = mprf__load_model_predictions_gui_devel_02;
catch
    
end
end








