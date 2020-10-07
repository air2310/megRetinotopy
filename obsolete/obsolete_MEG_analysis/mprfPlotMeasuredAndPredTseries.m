function mprfPlotMeasuredAndPredTseries(data, pred, periods, channels, gain)

if ~exist('gain','var') || isempty(gain)
    have_gain = false;
else
    have_gain = true;
    
end

cur_periods = 1:size(data,1);

if length(size(data)) == 3
    plot_data = nanmean(data,3);
    plot_data_ste = nanstd(data,[],3)./ sqrt(size(data,3));

else
    plot_data = data;
    plot_data_ste = zeros(size(data));
    
    
end

if length(size(pred)) == 3 
    plot_pred = nanmean(pred,3);
    plot_pred_ste = nanstd(pred,[],3)./ sqrt(size(pred,3));

else
    plot_pred = pred;
    plot_pred_ste = zeros(size(pred));

end





for n = 1:length(channels);
    
    % TODO Rescale predictions and data first:
    
    %m_pred = max(plot_pred(:,channels(n)));
    %m_data = max(plot_data(:,channels(n)));
    %ratio = m_pred ./ m_data;
    
    cur_plot_pred = plot_pred(:,channels(n));
    cur_plot_data = plot_data(:,channels(n));
    
  
    
    sel = ~isnan(cur_plot_pred) & ~isnan(cur_plot_data);
    
    cur_plot_pred(sel) = zscore(cur_plot_pred(sel));
    cur_plot_data(sel) = zscore(cur_plot_data(sel));
    
    if have_gain
        
        cur_plot_gain = abs(gain.M_02(:,channels(n)));
        cur_plot_gain(sel) = zscore(cur_plot_gain(sel));
        
    else
        cur_plot_gain = nan(size(cur_plot_pred));
        
    end
    
    y_lims = [-2 3];
    y_lims_02 = [-2 3];
    
    f=  figure; hold on;
    set(f,'Units','Normalized','Position',[0 1 1 .5]);
    set(gca,'FontSize',18,'FontWeight','bold','FontName','helvetica')
    
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
    
    
    patch(p_stim_periods,stim_y(:,ones(1,size(p_stim_periods,2))),'g')
    patch(p_blink_periods,blink_y(:,ones(1,size(p_blink_periods,2))),'r')
    patch(p_blank_periods,blank_y(:,ones(1,size(p_blank_periods,2))),'y')
    
    
    %p1 = plot(cur_periods, cur_plot_pred,'-','Color',[.7 .7 .7],'LineWidth',5);
    p3 = plot(cur_periods, cur_plot_gain,'-','Color',[.7 .7 .7],'LineWidth',5);

    
    p2 = plot(cur_periods, cur_plot_data,'k-','LineWidth',2);
    %     errorbar(cur_periods,plot_data_ste(:,n), plot_data_ste(:,n),'k.')
    plot(cur_periods,zeros(1,size(cur_plot_data,1)),'k--')
    legend([p2 p3],'Measured','Gain')
    
    cur_lim = axis;
%     errorbar(cur_periods,cur_plot_pred, plot_data_ste(:,n),'k.')
    
    axis(cur_lim);
    ylim(y_lims);
    title(sprintf('Channel nr: %d',channels(n)));
    xlabel('Stimulus period')
    ylabel('Zscored median sl_a_m_p / median bl_a_m_p')
    
end
























end