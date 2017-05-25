function mprfSessionRenderpRFActivityOnBrainstormSurface

thr_type = 'roi';

load(fullfile(pwd,'mprfSESSION.mat'));
global mprfSESSION
global fh

load(mprfSESSION.source.bs_head_model);
[~, surf.name] = fileparts(bs_model.SurfaceFile);
surf.type = 'brainstorm';


[pred_resp, stim] = mprfSessionGeneratePRFResponsePredictions(surf,[],'fmri',thr_type);

sz = size(stim.full_im);

fh = figure;
colormap gray;
params.ui.stim_slider = uicontrol('Style','slider',...
    'Units','Normalized',...
    'min',1,'max',sz(3),'value',1,...
    'SliderStep',[1/(sz(3)-1) 10/(sz(3)-1)],...
    'Position',[0.05 .9 .2 .1],...
    'Callback',@updateThePlot);


bs_msh = mprfMeshFromBrainstorm(mprfSESSION.source.(surf.name));
bs_msh = meshVisualize(bs_msh);


params.msh = bs_msh;
params.pred = pred_resp;
params.stim = stim.full_im;
params.cmap = jet(256);

set(gcf,'UserData',params);

updateThePlot;


end


function updateThePlot(hobj,junk)
global fh
params = get(fh,'UserData');
v = ceil(get(params.ui.stim_slider,'Value'));

stim_im = params.stim(:,:,v);
pred = params.pred(v,:);

pause(1/10); % Prevent mrMesh from 'suffocating'
doThePlotting(pred, stim_im, params.msh);


end


function doThePlotting(pred, image, msh)
global fh
params = get(fh,'UserData');



drange = [prctile(pred, 5) prctile(pred,95)];
msh = mprfSessionColorMesh(msh, pred, params.cmap, drange, ~isnan(pred));


cla;
imagesc(image);
axis image
axis off

params.msh = msh;
set(fh,'UserData',params);










end



