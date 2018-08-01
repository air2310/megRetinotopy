function fh = mprfPlotHeadLayout(high_light_channels,add_text,hl_size, plot_all_dots)

if ~exist('high_light_channels','var')
    high_light_channels = [];
end

if ~exist('add_text','var') || isempty(add_text)
    add_text = true;
end

if ~exist('hl_size','var') || isempty(hl_size)
    hl_size = 10;
end

if ~exist('plot_all_dots','var') || isempty(plot_all_dots)
    plot_all_dots = true;
end

fh = figure; hold on;


load(which('meg160_example_hdr.mat'))
layout = ft_prepare_layout([],hdr);

xpos = layout.pos(1:157,1);
ypos = layout.pos(1:157,2);
outline = layout.outline;

cur_chan = high_light_channels;


set(0,'CurrentFigure',fh);
for n = 1:length(outline)
    plot(outline{n}(:,1), outline{n}(:,2), 'k-','LineWidth',2)
end

if plot_all_dots
    plot(xpos(1:157),ypos(1:157),'ko','MarkerFaceColor','k')
end
if ~isempty(cur_chan)
    plot(xpos(cur_chan),ypos(cur_chan),'rd',...
        'MarkerFaceColor','r',...
        'MarkerSize',hl_size)
end

if add_text
    text(xpos(1:157)',ypos(1:157)',num2str((1:157)'));
end

axis tight
axis square
axis off


end