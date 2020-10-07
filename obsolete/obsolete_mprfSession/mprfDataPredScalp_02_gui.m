function mprfDataPredScalp_02_gui(data, pred, abs_pred, stim)

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
params_02.abs_pred_im = abs_pred;


set(pd_scalp_im,'UserData',params_02);

update_the_plot;

end


function update_the_plot(hobj,junk)

global pd_scalp_im
params = get(pd_scalp_im,'UserData');


cur_idx = round(get(params.ui.slider,'Value'));
cur_im = params.data_im(:,:,:,cur_idx);
cur_im_pred = params.pred_im(:,:,:,cur_idx);
cur_im_abs_pred = params.abs_pred_im(:,:,:,cur_idx);
cur_stim = params.stim(:,:,cur_idx);


set(params.ui.slider_text,'String',num2str(cur_idx));

set(0,'CurrentFigure',pd_scalp_im);

stim_asp_ratio = size(cur_stim,1) ./ size(cur_stim,2) ;
im_asp_ratio = size(cur_im,1) ./ size(cur_im,2); 

stim_margin_02 = 0.31;
stim_margin_01 = 0.02;
stim_height = .2;
stim_width = stim_height ./ stim_asp_ratio;

im_left = 0.25;
im_bottom = .02;
im_width = .5;
im_height = im_width .* im_asp_ratio;

abs_pred_im_left = 0.02;
abs_pred_im_bottom = .4;
abs_pred_im_width = .45;
abs_pred_im_height = abs_pred_im_width.* im_asp_ratio;

pred_im_left = 0.52;
pred_im_bottom = .4;
pred_im_width = .45;
pred_im_height = pred_im_width .* im_asp_ratio;



stim_sp = subplot('Position',...
    [stim_margin_02 (1-stim_margin_01 - stim_height) stim_width stim_height]);
stim_ax = imagesc(cur_stim(:,:,ones(1,1,3)).*255);
axis off

im_sp = subplot('Position',...
    [im_left im_bottom im_width im_height]);
im_ax = imagesc(cur_im);
axis off

im_sp = subplot('Position',...
    [abs_pred_im_left abs_pred_im_bottom abs_pred_im_width abs_pred_im_height]);
im_ax = imagesc(cur_im_abs_pred);
axis off


im_sp = subplot('Position',...
    [pred_im_left pred_im_bottom pred_im_width pred_im_height]);
im_ax = imagesc(cur_im_pred);
axis off

set(pd_scalp_im,'UserData',params);





end





























