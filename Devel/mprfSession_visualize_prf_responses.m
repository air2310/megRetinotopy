


global mprfSESSION
if isempty(mprfSESSION)
    load('mprfSESSION.mat')
end

pred = mprf__load_model_predictions;

bs_msh = mprfMeshFromBrainstorm(mprfSESSION.source.(surf.name));
bs_msh = meshVisualize(bs_msh);

[~,tmp] = fileparts(pred.bs.model_file);
surf_file = [tmp(strfind(tmp,'tess'):end) '.mat'];






















