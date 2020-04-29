% coherence may not be the best measure: just divide by neighbouring
% channels. Also, assumes linearity of response (i.e. no saturation).
% Assumes that increases in broadband signal are not happening, i.e. if
% broadband increases, denominator changes as well at other frequencies.
% Better to have a cut off frequency 


sub_004 = load('current_wlsubj_004/Plot_data.mat');
sub_040 = load('current_wlsubj_040/Plot_data.mat');

plot_sub = '004';
data_field = 'sl_co';
pred_field = 'M_02';

eval(['cur_data = sub_' plot_sub '.data_pp.' data_field ';']);
eval(['cur_pred = sub_' plot_sub '.gain.' pred_field ';']);

cur_data_ste = [];
try
    
    eval(['cur_data_ste = sub_' plot_sub '.pred_pp.' data_field '_ste;']);

catch
    
    
end

if strcmpi(pred_field, 'm_02')
    cur_pred = abs(cur_pred);
end


v_cur_data = nanvar(cur_data);
v_cur_pred = nanvar(cur_pred);

[s_v_cur_data, s_idx_cur_data] = sort(v_cur_data,'descend');

[s_v_cur_pred, s_idx_cur_pred] = sort(v_cur_pred,'descend');

fprintf('High values in data: %d\n',s_idx_cur_data(1:10))
fprintf('\n')
fprintf('High values in prediction: %d\n',s_idx_cur_pred(1:10))

fh_01 = figure;
fh_01 = megPlotMap(v_cur_data, [], fh_01,jet(256));
hgsave(fh_01, ['current_wlsubj_' plot_sub '/ScalpImage_data_' data_field]);


fh_02 = figure;
fh_02 = megPlotMap(v_cur_pred, [], fh_02,jet(256));
hgsave(fh_02, ['current_wlsubj_' plot_sub '/ScalpImage_pred_' pred_field]);

cur_r = nan(1,size(cur_pred,2));
cur_p = nan(1,size(cur_pred,2));


for n = 1:size(cur_pred,2)

    
    sel = ~isnan(cur_data(:,n)) & ~isnan(cur_pred(:,n)); 
    
    [tmp1, tmp2] = corrcoef(cur_data(sel,n), cur_pred(sel,n));
    cur_r(n) = tmp1(1,2);
    cur_p(n) = tmp2(1,2);

    
    
    
end

cur_r(isnan(cur_r)) = 0;

[s_v_cur_r, s_idx_cur_r] = sort(cur_r,'descend');
fprintf('\n')
fprintf('High values in correlation: %d\n',s_idx_cur_r(1:10))


fh_03 = figure;
fh_03 = megPlotMap(cur_r, [-1 1], fh_03,jet(256));
hgsave(fh_03, ['current_wlsubj_' plot_sub '/ScalpImage_corr_' data_field '_' pred_field]);

mprfPlotHeadLayout([]);


periods.blank = [3:5 30:32 57:59 84:86 111:113 138:140];
periods.blink = [1 2 28 29 55 56 82 83 109 110 136 137];
periods.stim = setdiff(1:140,[periods.blink periods.blank]);


channels = s_idx_cur_r(1:5);
head_plot_fh = nan(size(channels));
tseries_f = nan(size(channels));
plot_pred = -1;


for n = 1:length(channels);
    
    min_y = 0;
    
    if plot_pred == 1
        % Just prediction
        cur_plot = cur_pred(:,channels(n));
        
        y_lims = [min_y nanmax(cur_plot)];
        y_lims_02 = y_lims;
        
    elseif plot_pred == 0
        % Just data
        cur_plot = cur_data(:,channels(n));
        
        if ~isempty(cur_data_ste)
            cur_plot_ste = cur_data_ste(:,channels(n));
        end
            
        
        y_lims = [min_y nanmax(cur_plot)];
        y_lims_02 = y_lims;
        
    elseif plot_pred == -1
        % Both data and prediction
        tmp_pred = cur_pred(:,channels(n));
        tmp_data = cur_data(:,channels(n));
        
        
        sel = ~isnan(tmp_pred);
        tmp_pred(sel) = zscore(tmp_pred(sel));
        
        sel = ~isnan(tmp_data);
        std_data = std(tmp_data(sel));
        tmp_data(sel) = zscore(tmp_data(sel));
        
        
        if ~isempty(cur_data_ste)
            cur_plot_ste = (cur_data_ste(:,channels(n))./std_data);
            
        end
        
        
        cur_plot = [tmp_pred tmp_data];
        
        max_y = ceil(max(cur_plot(:) .* .5) ./ .5);
        
        y_lims = [-2.5 max_y];
        y_lims_02 = [-2.5 max_y];
        
    end
    cur_periods = 1:size(cur_plot,1);
    
    
    tseries_f(n) =  figure; hold on;
    set(tseries_f(n),'Units','Normalized','Position',[0 1 1 .5]);
    set(gca,'FontSize',24,'FontWeight','bold','FontName','helvetica')
    
    p_stim_periods = [reshape(periods.stim,[],5); ...
        flipud(reshape(periods.stim,[],5))];
    
    stim_y = y_lims(ones(size(p_stim_periods,1)/2,1),:);
    stim_y = stim_y(:);
    
    p_blink_periods = [reshape(periods.blink,[],6); ...
        flipud(reshape(periods.blink,[],6))];
    
    
    blink_y = y_lims(ones(size(p_blink_periods,1)/2,1),:);
    blink_y = blink_y(:);
    
    p_blank_periods = [reshape(periods.blank,[],6); ...
        flipud(reshape(periods.blank,[],6))];
    
    blank_y = y_lims(ones(size(p_blank_periods,1)/2,1),:);
    blank_y = blank_y(:);
    
    
    %patch(p_stim_periods,stim_y(:,ones(1,size(p_stim_periods,2))),'g')
    patch(p_blink_periods,blink_y(:,ones(1,size(p_blink_periods,2))),'r')
    patch(p_blank_periods,blank_y(:,ones(1,size(p_blank_periods,2))),'y')
    
    if size(cur_plot,2) == 2
        
        
        if ~isempty(cur_plot_ste)
            x_patch = [cur_periods fliplr(cur_periods)];
            y_patch = [(cur_plot(:,2) - cur_plot_ste)' flipud((cur_plot(:,2) + cur_plot_ste))'];
            sel = ~isnan(y_patch);
            p = patch(x_patch(sel), y_patch(sel),'b');
            set(p,'FaceColor',[.65 .65 .65])
            
        end
        
        
        pred_p = plot(cur_plot(:,1),'r','LineWidth',5);
        data_p = plot(cur_plot(:,2),'k-','LineWidth',4);
        
        legend([data_p pred_p],{'Measured','Predicted'})
        ylim(y_lims);
        
    elseif size(cur_plot,2) == 1 && ~plot_pred
        
        if ~isempty(cur_plot_ste)
            x_patch = [cur_periods fliplr(cur_periods)];
            y_patch = [(cur_plot - cur_plot_ste)' flipud((cur_plot + cur_plot_ste))'];
            sel = ~isnan(y_patch);
            p = patch(x_patch(sel), y_patch(sel),'b');
            set(p,'FaceColor',[.65 .65 .65])
            
        end
        plot(cur_periods, cur_plot, 'k-', 'LineWidth',4)
        
        ylim([min_y y_lims(2)]);
        
    elseif size(cur_plot,2) == 1 && plot_pred
        plot(cur_periods, cur_plot, 'color',[.5 .5 .5], 'LineWidth',4)
        ylim([min_y y_lims(2)]);
        
    end
    %p1 = plot(cur_periods, cur_plot_pred,'-','Color',[.7 .7 .7],'LineWidth',5);
    %p3 = plot(cur_periods, cur_plot_gain,'-','Color',[.7 .7 .7],'LineWidth',5);

    
    %p2 = plot(cur_periods, cur_plot_data,'k-','LineWidth',2);
    %     errorbar(cur_periods,plot_data_ste(:,n), plot_data_ste(:,n),'k.')
    %plot(cur_periods,zeros(1,size(cur_plot_data,1)),'k--')
    %legend([p2 p3],'Measured','Gain')
    
    cur_lim = axis;
%     errorbar(cur_periods,cur_plot_pred, plot_data_ste(:,n),'k.')
    
    axis(cur_lim);
    title(sprintf('Channel nr: %d',channels(n)));
    xlabel('Stimulus period')
    
    head_plot_fh(n) = mprfPlotHeadLayout(channels(n),false,20, false);
    set(head_plot_fh(n),'UserData',channels(n))
    
    cur_plot_ste = [];
end


% CO 004
%  0.5982    0.5773    0.5074    0.5029    0.4905    0.4555    0.4331    0.4245    0.4093    0.3907
% AMP 040
% 0.5863    0.5605    0.5005    0.4727    0.4537    0.4350    0.4288    0.4229    0.4209    0.3867
% C0 040
%  0.5926    0.4543    0.4339    0.4307    0.4303    0.4133    0.4053    0.4041    0.3985    0.3947
% AMP 040
%  0.5857    0.4330    0.4077    0.4077    0.3735    0.3725    0.3688    0.3674    0.3669    0.3651

return
do_figs = 1:5


for n = do_figs
    ppos = [100 100/3.5556];
    cur_fh = tseries_f(n);
    
    set(0,'CurrentFigure',tseries_f(n))
    
    set(gca,'color',[1 1 1],...
        'fontname','helvetica',...
        'fontweight','bold',...
        'fontsize',30,...
        'xticklabel','')
    
    xlabel('');
    
    set(cur_fh,'color','w',...
        'PaperUnits','centimeter',...
        'PaperSize',ppos,...
        'PaperPosition',[.1 1 ppos-[.2 2]])
    
    print(cur_fh, '-dpdf', ['tseries_fig_' num2str(cur_fh)])
    
    set(0,'CurrentFigure',head_plot_fh(n))
    set(head_plot_fh(n),'color','w')
    export_fig(['head_plot_fig_' num2str(cur_fh) '.png'],...
        '-transparent',true)
    
    
end


[s_v_cur_r(1:5);s_idx_cur_r(1:5)]


return

cur_fh = figure; 
cur_tseries_pred = cur_pred(:,channels(1));
pred_std = std(cur_pred(:,channels(1)))


cur_tseries_data = cur_pred(:,channels(1)) + randn(size(cur_tseries_pred)) * (pred_std /2);

plot(cur_tseries_data)
ppos = [100 100/3.5556];

set(0,'CurrentFigure',cur_fh)

set(gca,'color',[1 1 1],...
    'fontname','helvetica',...
    'fontweight','bold')

set(get(gca,'children'), ...
    'color','k',...
    'linestyle','-',...
    'LineWidth',6) 

axis off
axis tight

set(cur_fh,'color','w',...
    'PaperUnits','centimeter',...
    'PaperSize',ppos,...
    'PaperPosition',[.1 1 ppos-[.2 2]])

print -dpdf cur_fh


%    'color',[.5 .5 .5], ...
return

[stim_file, stim_path] = uigetfile('*','Please select stimulus used to generate the predictions/acquire the MEG data');
load(fullfile(stim_path, stim_file))

cur_pred = sub_040.gain.M_02;
cur_data = sub_040.data_pp.sl_co;
cur_stim = pred_stim.full_im;
% megPlotMap(sensor_data, clims, fh, cm , ttl, data_hdr,cfg)

cur_pred(isnan(cur_pred)) = 0;
cur_data(isnan(cur_data)) = 0;

pred_lim = max(abs(cur_pred(:)));
pred_lim = [-pred_lim pred_lim];

pred_lim_02 = [0 max(abs(cur_pred(:)))];


fh_pred = figure;
fh_data = figure;
fh_pred_02 = figure;

for n = 1:size(cur_pred,1)

    set(0,'CurrentFigure',fh_pred);
    fh_pred = megPlotMap(cur_pred(n,:),pred_lim,fh_pred,jet(256), sprintf('Prediction epoch %d',n))
    drawnow;
    tmp  = getframe(fh_pred);
    if n == 1
        pred_im = zeros([size(tmp.cdata) size(cur_pred,1)],'uint8');
        
    end
    pred_im(:,:,:,n) = tmp.cdata;
    
        
    set(0,'CurrentFigure',fh_pred_02);
    fh_pred_02 = megPlotMap(abs(cur_pred(n,:)),pred_lim_02,fh_pred_02,jet(256), sprintf('Absolute prediction epoch %d',n))
    drawnow;
    tmp  = getframe(fh_pred_02);
    if n == 1
        pred_im_02 = zeros([size(tmp.cdata) size(cur_pred,1)],'uint8');
        
    end
    pred_im_02(:,:,:,n) = tmp.cdata;
    
    
    set(0,'CurrentFigure',fh_data);
    fh_data = megPlotMap(cur_data(n,:),[0 .5],fh_data,jet(256), sprintf('Data epoch %d',n))
    drawnow;
    tmp  = getframe(fh_data);
    if n == 1
        data_im = zeros([size(tmp.cdata) size(cur_data,1)],'uint8');
        
    end
    data_im(:,:,:,n) = tmp.cdata;
    
    


end
