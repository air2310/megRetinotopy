function mprfDataPredScalp_gui(data, pred, stim)

global pd_scalp_im
pd_scalp_im = figure;

start_val = 1;
params_02.ui.slider = uicontrol('style','slider',...
    'Units','normalized',...
    'Position',[.05 .925 .2 .05],...
    'min',0,'max',size(stim,3),'value',start_val,...
    'SliderStep',[1/(size(stim,3)-1) 10/(size(stim,3)-1)],...
    'Callback',@update_the_plot);


params_02.ui.slider_text = uicontrol('style','text',...
    'Units','normalized',...
    'Position',[0.27 0.925, 0.05 0.05],...
    'string',num2str(start_val));
    


params_02.stim = uint8(stim./max(stim(:)));
params_02.data_im = data;
params_02.pred_im = pred;


set(pd_scalp_im,'UserData',params_02);

update_the_plot;

end


function update_the_plot(hobj,junk)

global pd_scalp_im
params = get(pd_scalp_im,'UserData');


cur_idx = round(get(params.ui.slider,'Value'));
cur_im = params.data_im(:,:,:,cur_idx);
cur_im_pred = params.pred_im(:,:,:,cur_idx);
cur_stim = params.stim(:,:,cur_idx);


set(params.ui.slider_text,'String',num2str(cur_idx));

set(0,'CurrentFigure',pd_scalp_im);

stim_asp_ratio = size(cur_stim,1) ./ size(cur_stim,2) ;
im_asp_ratio = size(cur_im,1) ./ size(cur_im,2); 

stim_margin_02 = 0.31;
stim_margin_01 = 0.02;
stim_height = .2;
stim_width = stim_height ./ stim_asp_ratio;

im_margin_01 = 0.02;
im_margin_02 = 0.02;
im_heigth = ((1-stim_margin_01 - stim_height) - im_margin_01*2) ./2;
im_width = (im_heigth ./ im_asp_ratio);

pred_im_margin_01 = 0.02;
pred_im_margin_02 = im_width;
pred_im_heigth = ((1-pred_im_margin_01 - stim_height) - pred_im_margin_01*2) ./2;
pred_im_width = (pred_im_heigth ./ im_asp_ratio);


stim_sp = subplot('Position',...
    [stim_margin_02 (1-stim_margin_01 - stim_height) stim_width stim_height]);
stim_ax = imagesc(cur_stim(:,:,ones(1,1,3)).*255);
axis off

im_sp = subplot('Position',...
    [im_margin_02 .25 im_width im_heigth]);
im_ax = imagesc(cur_im);
axis off

im_sp = subplot('Position',...
    [im_margin_02+pred_im_margin_02 .25 pred_im_width pred_im_heigth]);
im_ax = imagesc(cur_im_pred);
xlabel('Prediction')
axis off




set(pd_scalp_im,'UserData',params);





end





























