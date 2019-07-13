function mprfSessionRenderMEGPredictionOnScalpSurface

load(fullfile(pwd,'mprfSESSION.mat'));

if exist('ft_defaults','file')
else
    tbUse('retMEG')
end


[stim_im,head_im,M_02,M]  = mprfSessionGetPredictedData;

f1 = figure;
megPlotMap(mean(abs(M_02),1),[],f1,jet(256));
mprfSessionRenderMEGPredictionOnScalpSurface_gui(head_im, stim_im)



end

function mprfSessionRenderMEGPredictionOnScalpSurface_gui(head_im, stim_im)

global mprf_fh_02
mprf_fh_02 = figure;

start_val = 1;
params_02.ui.slider = uicontrol('style','slider',...
    'Units','normalized',...
    'Position',[.05 .925 .2 .05],...
    'min',0,'max',size(stim_im,3),'value',start_val,...
    'SliderStep',[1/(size(stim_im,3)-1) 10/(size(stim_im,3)-1)],...
    'Callback',@update_the_plot);


params_02.ui.slider_text = uicontrol('style','text',...
    'Units','normalized',...
    'Position',[0.27 0.925, 0.05 0.05],...
    'string',num2str(start_val));
    


params_02.stim = uint8(stim_im ./ max(stim_im(:)));
params_02.images = head_im;


set(mprf_fh_02,'UserData',params_02);

update_the_plot;

end


function update_the_plot(hobj,junk)

global mprf_fh_02
params = get(mprf_fh_02,'UserData');


cur_idx = round(get(params.ui.slider,'Value'));
cur_im = params.images(:,:,:,cur_idx);
cur_stim = params.stim(:,:,cur_idx);


set(params.ui.slider_text,'String',num2str(cur_idx));

set(0,'CurrentFigure',mprf_fh_02);

stim_asp_ratio = size(cur_stim,1) ./ size(cur_stim,2) ;
im_asp_ratio = size(cur_im,1) ./ size(cur_im,2); 

stim_margin_02 = 0.31;
stim_margin_01 = 0.02;
stim_height = .2;
stim_width = stim_height ./ stim_asp_ratio;

im_margin_01 = 0.02;
im_margin_02 = 0.02;
im_heigth = (1-stim_margin_01 - stim_height) - im_margin_01*2;
im_width = im_heigth ./ im_asp_ratio;

stim_sp = subplot('Position',...
    [stim_margin_02 (1-stim_margin_01 - stim_height) stim_width stim_height]);
stim_ax = imagesc(cur_stim(:,:,ones(1,1,3)).*255);
axis off


im_sp = subplot('Position',...
    [im_margin_02 im_margin_01 im_width im_heigth]);
im_ax = imagesc(cur_im);
axis off





set(mprf_fh_02,'UserData',params);





end












