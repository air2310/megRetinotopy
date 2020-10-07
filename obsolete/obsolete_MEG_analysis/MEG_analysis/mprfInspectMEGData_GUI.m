function mprfInspectMEGData_GUI

global mprf_fh
global data
% Blank periods
% Blink periods
% stim periods
% Needed:
% 1. be able to select the amount of channels to display and which
% channels. Edit field, with separate channles space separated
% 2. An option to compute the power spectrum of the presented data
% 3. Be able to select which repeat(s) and which period(s) to display
% Repeats same as channels. For the periods, I would like to have the
% stimulus next to it. Makes selection easier. Make a second figure? Like a
% montage?
% 4. Plot the current channels on a Scalp surface

if ~exist('data','var') || isempty(data);
    [data, source_file] = mprfGetGuiData;
    data.gui.source_file = source_file;
end

% at NYU or not:
nyu = true;

% Paths:
if nyu
    params.paths.data_dir = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wl_subj040/wl_subj040_20170406';
    params.paths.code_dir = '/Volumes/server/Projects/MEG/Retinotopy/Code';
    addpath(genpath('/Volumes/server/Projects/MEG/Retinotopy/Code'));
    
else
    params.paths.data_dir = '/Users/bpklein/Documents/Server_offline_folder/Data/MEG/wl_subj040/wl_subj040_20170406'; %#ok<UNRCH>
    params.paths.code_dir = '/Volumes/server/Projects/MEG/Retinotopy/Code';
    addpath(genpath('/Users/bpklein/Documents/Server_offline_folder/Code'));
end

params.paths.source_file = data.gui.source_file;
params.paths.raw_file = 'R1151_Retinotopy_04.06.17.sqd';                             % Raw file
params.paths.processed_file = 'Epoched_data_preproc.mat';                            % Epoched data file
params.paths.raw_mat_data_file = 'Raw_data_ft_read_data_data_channels.mat';          % Raw mat file, read with ft_read_data, data channels only
params.paths.raw_mat_tr_file = 'Raw_data_ft_read_data_trigger_channels.mat';         % Raw mat file, read with ft_read_data, trigger channels only
params.paths.raw_mat_ph_dile = 'Raw_data_ft_read_data_photo_diode_channels.mat';     % Raw mat file, read with ft_read_data, photodiode only
params.paths.param_dir = fullfile(params.paths.data_dir,'Stimulus/Param_files');                  % Parameter files
params.paths.stim_dir = fullfile(params.paths.data_dir,'Stimulus/Stimulus_files');
params.paths.example_hdr = fullfile(params.paths.code_dir,'mprfSession','external','meg160_example_hdr.mat');
params.paths.amp_head_pp_file = '';

good_channels.sweep_0102 = [25    51    69    13    26    43    88     2     1    60];

rnd_idx = randperm(length(good_channels.sweep_0102));
def_channels = num2str(good_channels.sweep_0102(rnd_idx(1:5)));

mprf_fh = figure;

box_d = 0.06;

%% Channels edit field
box_n = 1;
params.channels.ui_edit_01 = uicontrol(...
    'Style','edit',...
    'Units','Normalized',...
    'String',def_channels,...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@update_plot_params);

%% Channels text
params.channels.ui_text_01 = uicontrol(...
    'Style','text',...
    'Units','normalized',...
    'Position',[.05 (1 - box_n * box_d) .15 .05],...
    'String','Which channels?');

%% Repeats edit field
box_n = 2;
params.repeats.ui_edit_01 = uicontrol(...
    'Style','edit',...
    'Units','Normalized',...
    'String','1:5',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@update_plot_params);

%% Repeats text
params.repeats.ui_text_01 = uicontrol(...
    'Style','text',...
    'Units','normalized',...
    'Position',[.05 (1 - box_n * box_d) .15 .05],...
    'String','Which repeats?');

%% Repeats periods edit
box_n = 3;
params.periods.ui_edit_01 = uicontrol(...
    'Style','edit',...
    'Units','Normalized',...
    'String','6:27',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@update_plot_params);

%% Periods text:
params.periods.ui_text_01 = uicontrol(...
    'Style','text',...
    'Units','normalized',...
    'Position',[.05 (1 - box_n * box_d) .15 .05],...
    'String','Which periods?');


%% Stimulus figure pushbutton
box_n = 4;
params.stim.ui_pb_01 = uicontrol(...
    'Style','pushbutton',...
    'Units','Normalized',...
    'String','Stimulus figure',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@make_stimulus_figure);


%% Time series plot
box_n = 5;
params.plot.pb_make_plot = uicontrol(...
    'Style','pushbutton',...
    'Units','Normalized',...
    'String','Make plot',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@do_the_plotting);


%% PS current plot pushbutton
box_n = 6;
params.plot.pb_make_plot = uicontrol(...
    'Style','pushbutton',...
    'Units','Normalized',...
    'String','Plot powerspectrum of figure',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@power_spect_current_figure);

%% Head layout current channels:
box_n = 7;
params.plot.pb_make_ps_plot = uicontrol(...
    'Style','pushbutton',...
    'Units','Normalized',...
    'String','Make head layout',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@make_head_plot);


%% Plot average of figure
box_n = 8;
params.plot.pb_make_av_plot = uicontrol(...
    'Style','pushbutton',...
    'Units','Normalized',...
    'String','Plot average of figure',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@make_average_of_figure);


%% Histogram of SNR
box_n = 9;
params.plot.pb_make_snr_histogram = uicontrol(...
    'Style','pushbutton',...
    'Units','Normalized',...
    'String','Make SNR histogram',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@make_snr_histogram);

%% Amplitude at stimulus frequency across head
box_n = 10;
params.plot.make_amp_contrasts = uicontrol(...
    'Style','pushbutton',...
    'Units','Normalized',...
    'String','Amplitude contrasts',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@make_amp_contrasts);

%% Amplitude at stimulus frequency across head per periods
box_n = 11;
params.plot.make_amp_head_layout_pp = uicontrol(...
    'Style','pushbutton',...
    'Units','Normalized',...
    'String','Amplitude across head per position',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@make_amp_head_layout_pp);

%% Align phase
box_n = 12;
params.plot.align_phase = uicontrol(...
    'Style','pushbutton',...
    'Units','Normalized',...
    'String','Align phase of figure data',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@align_phase);

box_n = 13;
params.plot.plot_tseries = uicontrol(...
    'Style','pushbutton',...
    'Units','Normalized',...
    'String','Plot t-series',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@plot_tseries);

%% Amplitude at stimulus frequency across head per periods
box_n = 14;
params.plot.plot_amp_tseries = uicontrol(...
    'Style','pushbutton',...
    'Units','Normalized',...
    'String','Plot amplitude t-series',...
    'Position',[0.26 (1 - box_n * box_d) .5 .05],...
    'Callback',@plot_amp_tseries);

%% Store data:
params.periods.blank = [3:5 30:32 57:59 84:86 111:113 138:140];
params.periods.blink = [1 2 28 29 55 56 82 83 109 110 136 137];
params.periods.stim = setdiff(1:140,[params.periods.blink params.periods.blank]);



params.data = data;
set(mprf_fh,'UserData',params)

%% Set the default values in the parameter struct:
update_plot_params(params.channels.ui_edit_01);
update_plot_params(params.periods.ui_edit_01);
update_plot_params(params.repeats.ui_edit_01);



end

function make_stimulus_figure(hobj,junk) %#ok<*INUSD>
%% To create the stimulus figure

global mprf_fh
params = get(mprf_fh,'UserData');


create_stim_figure = true;
if isfield(params,'stim')
    if isfield(params.stim,'stim_figure')
        if exist([params.stim.stim_figure '.fig'],'file')
            try
                stim_fh = hgload(params.stim.stim_figure);
                create_stim_figure = false;
                
                
            catch
                
            end
        end
    end
end

if create_stim_figure
    
    stimfiles = dir(fullfile(params.paths.stim_dir,'*.mat'));
    stim = load(fullfile(params.paths.stim_dir, stimfiles(1).name));
    
    stim_range = double([min(stim.stimulus.images(:)) max(stim.stimulus.images(:))]);
    bk = double(mode(mode(stim.stimulus.images(:,:,1))));
    new_period = stim.stimulus.seq(stim.stimulus.trigSeq > 0);
    
    im_out_size = [101 101];
    im_seq = zeros([im_out_size, length(new_period)],'uint8');
    
    
    for n = 1:length(new_period)
        tmp_im = double(stim.stimulus.images(:,:,new_period(n)));
        tmp_im = imresize(tmp_im,im_out_size,'nearest');
        im_seq(:,:,n) = uint8(ceil(abs(tmp_im - bk) ./ max(abs(stim_range - bk))));
        
    end
    
    if length(new_period) == 140
        nrows = 10;
        ncols = 14;
    end
    
    stim_fh = figure;
    colormap gray
    
    for n = 1:nrows
        for nn = 1:ncols
            pos = ((n-1) * ncols) + nn;
            
            s = subplot(nrows,ncols,pos);
            set(s,'pos',get(s,'pos').*[1 1 1.05 1.05]);
            imagesc(im_seq(:,:,pos),[0 1])
            title(sprintf('%d',pos))
            xlabel(pos)
            axis off
            axis image
        end
    end
    
    params.stim.stim_figure = tempname;
    hgsave(stim_fh, params.stim.stim_figure);
    
    
end

params.stim.stim_figure_fh = stim_fh;
set(mprf_fh, 'UserData',params);

end

function update_plot_params(hobj,junk)
%% Set the channels and/or repeats:

global mprf_fh
params = get(mprf_fh,'UserData');


if hobj == params.channels.ui_edit_01
    val = get(params.channels.ui_edit_01,'String');
    cel_val = strsplit(val,' ');
    params.plot.channels = unique(cell2mat(cellfun(@eval,cel_val,'UniformOutput',false)));
    
elseif hobj == params.periods.ui_edit_01
    val = get(params.periods.ui_edit_01,'String');
    if strcmpi(val,'stimulus')
        params.plot.periods = params.periods.stim;
    elseif strcmpi(val,'blink')
        params.plot.periods = params.periods.blink;
    elseif strcmpi(val,'blank')
        params.plot.periods = params.periods.blank;
        
    else
        cel_val = strsplit(val,' ');
        params.plot.periods = unique(cell2mat(cellfun(@eval,cel_val,'UniformOutput',false)));
        
        
        
    end
    
elseif hobj == params.repeats.ui_edit_01
    val = get(params.repeats.ui_edit_01,'String');
    cel_val = strsplit(val,' ');
    params.plot.repeats = unique(cell2mat(cellfun(@eval,cel_val,'UniformOutput',false)));
    
    
else
    
    
end


set(mprf_fh,'UserData',params);



end

%% Plot time series

function do_the_plotting(hobj,junk)

global mprf_fh
global data
params = get(mprf_fh,'UserData');

cur_data = data.data(:,params.plot.periods, params.plot.repeats,params.plot.channels);

% plot different periods as a continous line (2), plot different repeats
% together (3) with separate plots for every channel (4) (subplot)

text_xpos = (size(cur_data,1) .* (0:size(cur_data,2)-1)) + size(cur_data,1)/2;



if size(cur_data,2) > 1
    cur_data = reshape(cur_data,[size(cur_data,1) .* size(cur_data,2), size(cur_data,3), size(cur_data,4)]);
else
    cur_data = squeeze(cur_data);
end

fh = figure;

for n = 1:size(cur_data,3)
    subplot(size(cur_data,3),1,n);
    plot(cur_data(:,:,n));
    title(sprintf('Channel nr: %d',params.plot.channels(n)));
    ypos = max(max(cur_data(:,:,n)));
    text(text_xpos',ypos(:,ones(1,size(text_xpos,2)))',num2str(params.plot.periods'))
    
    if n == size(cur_data,3)
        xlabel('Samples')
        
    end
    
    ylabel('Tesla')
    
end

plot_data.cur_data = cur_data;
plot_data.channels = params.plot.channels;
plot_data.periods = params.plot.periods;
plot_data.repeats = params.plot.repeats;

set(fh, 'UserData',plot_data);
set(mprf_fh,'UserData',params);


end


%% Plot powerspectrum

function power_spect_current_figure(hobj, junk)



fh = figure;
if isobject(fh)
    cur_fh = str2double(inputdlg('Which figure?','',1,{num2str(fh.Number-1)}));
    
else
    cur_fh = str2double(inputdlg('Which figure?','',1,{num2str(fh-1)}));
    
end

plot_data = get(cur_fh,'UserData');
max_fs = 65;

set(gca,'FontSize',12);
plot(ones(2,size(plot_data.cur_data,2)))
c_order = get(gca,'ColorOrder');
cla;

bl_range = 5;
sl_range = .05;
stim_freq = 10;
bl_start_end = stim_freq + [-bl_range/2 bl_range/2];
sl_start_end = stim_freq + [-sl_range/2 sl_range/2];

freq_range = [bl_start_end(1) sl_start_end bl_start_end(2)];

fprintf('Baseline range: %0.2f : %0.2f\n',bl_start_end)
fprintf('Stimulus range: %0.2f : %0.2f\n',sl_start_end)

% av_amp = cell(1,size(plot_data.cur_data,3));

for n = 1:size(plot_data.cur_data, 2)
    for nn = 1:size(plot_data.cur_data, 3)
        
        
        tmp_data = plot_data.cur_data(~isnan(plot_data.cur_data(:,n,nn)),n,nn);
        if any(tmp_data)
            tmp = fft(tmp_data);
            tmp = tmp(1:floor(size(tmp,1)/2)+1,:,:);
            tmp = (2 .* abs(tmp) ./ size(tmp_data,1));
            
            
            end_idx = ceil(mprfFreq2Index(size(tmp_data,1),max_fs,1000));
            plot_x = linspace(0,max_fs,end_idx);
            
            idx_range = round(mprfFreq2Index(size(tmp_data,1),freq_range,1000));
            
            
            nb_idx = [idx_range(1):idx_range(2) idx_range(3) : idx_range(4)];
            
            sl_idx = idx_range(2) : idx_range(3);
            
            spec_amp_nb = abs(tmp(nb_idx));
            spec_amp_sl = abs(tmp(sl_idx));
            
            av_nb = exp(nanmean(log(spec_amp_nb)));
            av_sl = exp(nanmean(log(spec_amp_sl)));
            
            snr = av_sl ./ av_nb;
            fprintf('Amplitude modulation for channel %d: %0.2f\n',...
                plot_data.channels(nn), snr)
            
            
            subplot(size(plot_data.cur_data,3),1,nn)
            hold on;
            
            cl_idx = mod(n,size(c_order,1));
            if cl_idx == 0
                cl_idx = size(c_order,1);
            end
            
            
            plot(plot_x, log(tmp(1:end_idx).^2),'--','Color',c_order(cl_idx,:));
            cur_y_lim = ylim;
            
            bl_patch_x = plot_x([(idx_range(1) : idx_range(4)) (idx_range(4) :-1: idx_range(1))]);
            sl_patch_x = plot_x([sl_idx fliplr(sl_idx)]);
            
            bl_patch_y = cur_y_lim(ones(length(bl_patch_x)/2,1),:);
            bl_patch_y = bl_patch_y(:);
            sl_patch_y = cur_y_lim(ones(length(sl_patch_x)/2,1),:);
            sl_patch_y = sl_patch_y(:);
            
        end
        
        
        if n == size(plot_data.cur_data,2)
            title(sprintf('Channel nr: %d',plot_data.channels(nn)))
            sl_patch = patch(sl_patch_x, sl_patch_y,'r');
            bl_patch = patch(bl_patch_x, bl_patch_y,'k');
            set(bl_patch,'FaceColor',[.5 .5 .5])
            
            set(gca,'Children',flipud(get(gca,'Children')));
            
            % Do patch here
            
        end
        
        ylabel('Log power')
        
        if nn == size(plot_data.cur_data,3) && n == 1
            xlabel('Frequency')
        end
        
    end
    
end

plot_data.cur_data = [];
set(fh,'UserData',plot_data)



end


%% Make head plot

function make_head_plot(hobj,junk)



fh = figure; hold on;


global mprf_fh
params = get(mprf_fh,'UserData');

get_layout_data = true;
if isfield(params,'plot')
    if isfield(params.plot,'head')
        if isfield(params.plot.head,'xpos')
            get_layout_data = false;
            
        end
    end
end

if get_layout_data
    load(params.paths.example_hdr)
    layout = ft_prepare_layout(hdr);
    
    params.plot.head.xpos = layout.pos(1:157,1);
    params.plot.head.ypos = layout.pos(1:157,2);
    params.plot.head.outline = layout.outline;
end


cur_chan = params.plot.channels;


set(0,'CurrentFigure',fh);
for n = 1:length(params.plot.head.outline)
    plot(params.plot.head.outline{n}(:,1), params.plot.head.outline{n}(:,2), 'k-','LineWidth',2)
end

plot(params.plot.head.xpos(1:157),params.plot.head.ypos(1:157),'ko','MarkerFaceColor','k')
plot(params.plot.head.xpos(cur_chan),params.plot.head.ypos(cur_chan),'rd','MarkerFaceColor','r')
text(params.plot.head.xpos(1:157)',params.plot.head.ypos(1:157)',num2str((1:157)'))
axis square
axis off

set(mprf_fh,'UserData',params);


end

%% Make average figure

function make_average_of_figure(hobj, junk)

global mprf_fh
params = get(mprf_fh,'UserData');

fh = figure; hold on;
if isobject(fh)
    cur_fh = str2double(inputdlg('Which figure?','',1,{num2str(fh.Number-1)}));
    
else
    cur_fh = str2double(inputdlg('Which figure?','',1,{num2str(fh-1)}));
    
end

plot_data = get(cur_fh,'UserData');

mean_signal = squeeze(nanmean(plot_data.cur_data,2));
std_signal = squeeze(nanstd(plot_data.cur_data,[],2));
ste_signal = std_signal ./ sqrt(size(plot_data.cur_data,2));


n_periods = size(plot_data.periods,2);
text_xpos = size(plot_data.cur_data,1) .* ((0:n_periods-1) ./ n_periods) ...
    + (size(plot_data.cur_data,1) / (2*n_periods));

for n = 1:size(mean_signal,2)
    subplot(size(mean_signal,2),1,n);
    hold on;
    
    cur_x = find(std_signal(:,n)~=0 & ~isnan(std_signal(:,n)))';
    cur_y = [mean_signal(cur_x,n) - std_signal(cur_x,n); flipud(mean_signal(cur_x,n) + std_signal(cur_x,n))];
    cur_x = [cur_x';flipud(cur_x')];
    
    %     tmp_p = patch(cur_x,cur_y,'k');
    %     set(tmp_p,'FaceColor',[.5 .5 .5],'EdgeColor','none');
    plot(mean_signal(:,n),'k-','LineWidth',2);
    
    
    
    title(sprintf('Channel nr: %d',params.plot.channels(n)));
    ypos = max(max(std_signal(:,n)));
    text(text_xpos',ypos(:,ones(1,size(text_xpos,2)))',num2str(plot_data.periods'))
    
    if n == size(plot_data.cur_data,2)
        xlabel('Samples')
        
    end
    
    ylabel('Tesla')
    
end


plot_data.cur_data = [];
plot_data.cur_data(:,1,:) = mean_signal;
set(fh,'UserData',plot_data);



end

%% Make phase plot
function make_snr_histogram(hobj, junk)

global mprf_fh
global data


params = get(mprf_fh, 'Userdata');

cur_periods = sort([params.periods.blank params.periods.stim]);

stim_freq = 10;
bl_freq = [9 11];

nb_idx_tmp = round(mprfFreq2Index(size(data.data,1),bl_freq,1000));
sl_idx_tmp = round(mprfFreq2Index(size(data.data,1),stim_freq,1000));

blank_snr = nan(numel(params.periods.blank), size(data.data,3) ,size(data.data,4));
stim_snr = nan(numel(params.periods.stim), size(data.data,3) ,size(data.data,4));

for nn = 1:size(data.data, 4);
    
    cur_blank_idx = 0;
    cur_stim_idx = 0;
    
    for n = 1:length(cur_periods)
        
        
        
        cur_data = squeeze(data.data(:,cur_periods(n),: , nn));
        
        if any(~isnan(cur_data(:)))
            
            ft_data = fft(cur_data);
            amp_data = (2 .* abs(ft_data) ./ size(ft_data,1));
            snr = exp(nanmean(log(amp_data(sl_idx_tmp,:)))) ...
                ./ exp(nanmean(log(amp_data(nb_idx_tmp,:))));
            
            
            if ismember(cur_periods(n), params.periods.blank)
                cur_blank_idx = cur_blank_idx + 1;
                blank_snr(cur_blank_idx, :, nn) = snr;
                
            elseif ismember(cur_periods(n), params.periods.stim)
                cur_stim_idx = cur_stim_idx + 1;
                stim_snr(cur_stim_idx, :, nn) = snr;
                
            else
                warning('Unknown periods');
                
            end
            
        end
        
        
        
        
    end
end

blank_snr_res = reshape(blank_snr,[],size(blank_snr,3));
stim_snr_res = reshape(stim_snr,[],size(stim_snr,3));

blank_snr_med = nanmedian(blank_snr_res,1);
stim_snr_med = nanmedian(stim_snr_res,1);

all_snr = stim_snr_med ./ blank_snr_med;
[all_snr_s, s_idx] = sort(all_snr, 'descend');

s_idx = s_idx(~isnan(all_snr_s));
all_snr_s = all_snr_s(~isnan(all_snr_s));

figure;
plot(s_idx, all_snr_s,'k-')

n_channels_to_plot = 5;
for n = 1:n_channels_to_plot;
    cur_chan = s_idx(n);
    
    [blank_y, blank_x] = hist(blank_snr_res(:,cur_chan), 0:.5:10.5);
    [stim_y, stim_x] = hist(stim_snr_res(:,cur_chan), 0:.5:10.5);
    
    blank_y = blank_y ./ sum(blank_y);
    stim_y = stim_y ./ sum(stim_y);
    
    figure;
    b_h = bar([blank_x' stim_x'],[blank_y', stim_y'],1);
    set(b_h(1), 'FaceColor', [0 0 0]);
    set(b_h(2), 'FaceColor', [1 1 1]);
    legend(b_h, 'Blank periods','Stimulus periods');
    ylabel('Proportion of periods')
    xlabel('SNR')
    title(sprintf('Channel nr %d',cur_chan));
    
end

fprintf('Best channels: \n');
fprintf('%d\n',s_idx(1:10)');

fprintf('\nWorst channels: \n');
fprintf('%d\n',s_idx(end-9:end)');



fh = figure;
megPlotMap(all_snr,[ceil(min(all_snr(:)./0.05)) .* 0.05 ceil(max(all_snr(:)./0.05)) .* 0.05],...
    fh,jet(256));


end


function make_amp_contrasts(hobj, junk)


global data
global mprf_fh

params = get(mprf_fh,'UserData');

stimfiles = dir(fullfile(params.paths.stim_dir,'*.mat'));
stim = load(fullfile(params.paths.stim_dir, stimfiles(1).name));

clear stimfiles

stim_range = double([min(stim.stimulus.images(:)) max(stim.stimulus.images(:))]);
bk = double(mode(mode(stim.stimulus.images(:,:,1))));
new_period = stim.stimulus.seq(stim.stimulus.trigSeq > 0);

im_out_size = [101 101];
im_seq = zeros([im_out_size, length(new_period)],'uint8');


for n = 1:length(new_period)
    tmp_im = double(stim.stimulus.images(:,:,new_period(n)));
    tmp_im = imresize(tmp_im,im_out_size,'nearest');
    im_seq(:,:,n) = uint8(ceil(abs(tmp_im - bk) ./ max(abs(stim_range - bk))));
    
end

blank_periods = params.periods.blank;
blink_periods = params.periods.blink;
stim_periods = params.periods.stim;

up_periods = [33:36 73:81 106:108 114:123];
low_periods = [66:68 87:90 60:67 129:135];
lhs_periods = [6:14 66:68 106:108];
rhs_periods = [19:26 33:36 87:90];

cur_ratio = numel(stim_periods) / numel(blank_periods);
stim_vs_blank = zeros(1,size(im_seq,3));
stim_vs_blank(stim_periods) = 1;
stim_vs_blank(blank_periods) = -1 .* cur_ratio;

cur_ratio = numel(lhs_periods) / numel(rhs_periods);
lhs_vs_rhs = zeros(1,size(im_seq,3));
lhs_vs_rhs(lhs_periods) = 1;
lhs_vs_rhs(rhs_periods) = -1 .* cur_ratio;

cur_ratio = numel(up_periods) / numel(low_periods);
up_vs_low = zeros(1,size(im_seq,3));
up_vs_low(up_periods) = 1;
up_vs_low(low_periods) = -1 .* cur_ratio;

cur_ratio = numel(up_periods) / numel(blank_periods);
up_vs_blank = zeros(1,size(im_seq,3));
up_vs_blank(up_periods) = 1;
up_vs_blank(blank_periods) = -1 .* cur_ratio;

cur_ratio = numel(low_periods) / numel(blank_periods);
low_vs_blank = zeros(1,size(im_seq,3));
low_vs_blank(low_periods) = 1;
low_vs_blank(blank_periods) = -1  .* cur_ratio;

cur_ratio = numel(lhs_periods) / numel(blank_periods);
left_vs_blank = zeros(1,size(im_seq,3));
left_vs_blank(lhs_periods) = 1;
left_vs_blank(blank_periods) = -1 .* cur_ratio;

cur_ratio = numel(rhs_periods) / numel(blank_periods);
right_vs_blank = zeros(1,size(im_seq,3));
right_vs_blank(rhs_periods) = 1;
right_vs_blank(blank_periods) = -1 .* cur_ratio;


stim_freq = 10;
bl_freq = [9 11];

fprintf('Computing amplitude modulation for %d channels:\n', size(data.data,4))

snr_temp = ones(size(data.data,2),size(data.data,4));

for n = 1:size(data.data,4)
    fprintf('%d.',n);
    
    tmp_data = data.data(:,:,:,n);
    
    if any(~isnan(tmp_data(:)))
        tmp = fft(tmp_data,[],1);
        tmp = tmp(1:floor(size(tmp,1)/2)+1,:,:);
        
        nb_idx = round(mprfFreq2Index(size(tmp_data,1),bl_freq,1000));
        sl_idx = round(mprfFreq2Index(size(tmp_data,1),stim_freq,1000));
        
        spec_amp_nb = abs(tmp(nb_idx,:,:));
        spec_amp_sl = abs(tmp(sl_idx,:,:));
        
        av_nb = squeeze(exp(nanmean(log(spec_amp_nb))));
        av_sl = squeeze(spec_amp_sl);
        
        snr_temp(:,n) = nanmedian(av_sl,2) ./ nanmedian(av_nb,2);
        
        
        
    end
    
    
    
    
    
end



fprintf('Done\n')
snr_temp(isnan(snr_temp)) = 0;


stim_vs_blank_temp = stim_vs_blank * snr_temp;
lhs_vs_rhs_temp = lhs_vs_rhs * snr_temp;
up_vs_low_temp = up_vs_low * snr_temp;

up_vs_blank_temp = up_vs_blank * snr_temp;
low_vs_blank_temp = low_vs_blank * snr_temp;
left_vs_blank_temp = left_vs_blank * snr_temp;
right_vs_blank_temp = right_vs_blank * snr_temp;


figure;
megPlotMap(stim_vs_blank_temp,[floor(min(stim_vs_blank_temp./0.05) .* 0.05) ceil(max(stim_vs_blank_temp./0.05) .* 0.05)],...
    [],jet(256),'Stim vs blank, temporal averaging');

figure;
megPlotMap(lhs_vs_rhs_temp,[floor(min(left_vs_blank_temp./0.05) .* 0.05) ceil(max(left_vs_blank_temp./0.05) .* 0.05)], ...
    [],jet(256),'left vs right, temporal averaging');

figure;
megPlotMap(up_vs_low_temp,[floor(min(up_vs_low_temp./0.05) .* 0.05) ceil(max(up_vs_low_temp./0.05) .* 0.05)], ...
    [],jet(256),'up vs low, temporal averaging');


figure;
megPlotMap(up_vs_blank_temp,[floor(min(up_vs_blank_temp./0.05) .* 0.05) ceil(max(up_vs_blank_temp./0.05) .* 0.05)],...
    [],jet(256),'Up vs blank, temporal averaging');


figure;
megPlotMap(low_vs_blank_temp,[floor(min(low_vs_blank_temp./0.05) .* 0.05) ceil(max(low_vs_blank_temp./0.05) .* 0.05)],...
    [],jet(256),'Low vs blank, temporal averaging');


figure;
megPlotMap(left_vs_blank_temp,[floor(min(left_vs_blank_temp./0.05) .* 0.05) ceil(max(left_vs_blank_temp./0.05) .* 0.05)], ...
    [],jet(256),'left vs blank, temporal averaging');

figure;
megPlotMap(right_vs_blank_temp,[floor(min(right_vs_blank_temp./0.05) .* 0.05) ceil(max(right_vs_blank_temp./0.05) .* 0.05)], ...
    [],jet(256),'right vs blank, temporal averaging');



end


function make_amp_head_layout_pp(hobj, junk)

global data
global mprf_fh

params = get(mprf_fh,'UserData');
source_dir = fileparts(params.paths.source_file);
cur_dir = pwd;
cd(source_dir);
make_gui_images = true;
if isfield(params,'paths')
    if isfield(params.paths,'amp_head_pp_file')
        
        if isempty(params.paths.amp_head_pp_file)
            [fname, fpath] = uigetfile('Please select file to plot');
            
            if fname
                load(fullfile(fpath,fname));
                make_gui_images = false;
                
                
            end
        else
            
            try
                load(params.paths.amp_head_pp_file)
                
                make_gui_images = false;
                
            catch
                
                fprintf('Could not load stored images, creating new ones\n');
                
            end
        end
        
        
    end
end

cd(cur_dir)

if make_gui_images
    stimfiles = dir(fullfile(params.paths.stim_dir,'*.mat'));
    stim = load(fullfile(params.paths.stim_dir, stimfiles(1).name));
    
    clear stimfiles
    
    stim_range = double([min(stim.stimulus.images(:)) max(stim.stimulus.images(:))]);
    bk = double(mode(mode(stim.stimulus.images(:,:,1))));
    new_period = stim.stimulus.seq(stim.stimulus.trigSeq > 0);
    
    im_out_size = [101 101];
    im_seq = zeros([im_out_size, length(new_period)],'uint8');
    
    
    for n = 1:length(new_period)
        tmp_im = double(stim.stimulus.images(:,:,new_period(n)));
        tmp_im = imresize(tmp_im,im_out_size,'nearest');
        im_seq(:,:,n) = uint8(ceil(abs(tmp_im - bk) ./ max(abs(stim_range - bk))));
        
    end
    
    cur_data = data.data(:,params.periods.stim,:,:);
    cur_data = reshape(cur_data,[],size(data.data,4));
    fprintf('Computing amplitude modulation for %d channels:\n', size(cur_data,2))
    
    stim_freq = 10;
    bl_freq = [9 11];
    
    snr = ones(size(data.data,2),size(data.data,4));
    %tmp_data = squeeze(nanmean(data.data,3));
    
    for n = 1:size(data.data,4)
        fprintf('%d.',n);
        
        tmp_data_02 = data.data(:,:,:,n);
        
        if any(~isnan(tmp_data_02(:)))
            tmp = fft(tmp_data_02,[],1);
            tmp = tmp(1:floor(size(tmp,1)/2)+1,:,:);
            
            nb_idx = round(mprfFreq2Index(size(tmp_data_02,1),bl_freq,1000));
            sl_idx = round(mprfFreq2Index(size(tmp_data_02,1),stim_freq,1000));
            
            spec_amp_nb = abs(tmp(nb_idx,:,:));
            spec_amp_sl = abs(tmp(sl_idx,:,:));
            
            av_nb = squeeze(exp(nanmean(log(spec_amp_nb))));
            av_sl = squeeze(spec_amp_sl);
            
            snr(:,n) = nanmedian(av_sl,2) ./ nanmedian(av_nb,2);
            %    fprintf('Amplitude modulation for channel %d: %0.2f\n',...
            %        plot_data.channels(nn), snr)
            
        end
        
    end
    
    
    
    fprintf('Done\n')
    
    
    fh_02 = figure;
    
    
    
    for n = 1:size(snr,1)
        
        fh_02 = megPlotMap(snr(n,:),[0.75 max(snr(:)./0.05) .* 0.05],fh_02,jet(256),sprintf('Stimulus position %d',n));
        drawnow;
        
        tmp  = getframe(fh_02);
        
        if n == 1
            head_im = zeros([size(tmp.cdata) size(snr,1)],'uint8');
            
        end
        head_im(:,:,:,n) = tmp.cdata;
        
        
    end
    
    cur_time = datestr(now);
    cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
    
    source_dir = fileparts(params.paths.source_file);
    
    save(fullfile(source_dir,['Amplitude_head_overlay_per_period_' cur_time]),'head_im','im_seq')
    
    
    params.paths.amp_head_pp_file = fullfile(params.paths.data_dir,['Amplitude_head_overlay_per_period_' cur_time]);
    
end




mprfSNRPerPeriodGUI(params, im_seq,head_im);


end



function align_phase(hobj, junk)

global data

fh = figure; hold on;
if isobject(fh)
    cur_fh = str2double(inputdlg('Which figure?','',1,{num2str(fh.Number-1)}));
    
else
    cur_fh = str2double(inputdlg('Which figure?','',1,{num2str(fh-1)}));
    
end

plot_data = get(cur_fh,'UserData');


cur_data = data.data(:,plot_data.periods,plot_data.repeats,plot_data.channels);

sz = [size(cur_data,1) size(cur_data,2) size(cur_data,3) size(cur_data,4)];



for nnn = 1:sz(4)
    break_out_of_loop = false;
    
    
    for n = 1:sz(2)
        for nn = 1:sz(3)
            
            tmp = cur_data(:,n,nn,nnn);
            if any(isnan(tmp))
                
            else
                F_idx = mprfFreq2Index(size(tmp,1),10,1000);
                
                ft_tmp = fft(tmp);
                cur_ph = angle(ft_tmp(F_idx));
                break_out_of_loop = true;
                
            end
            
            if break_out_of_loop
                break
            end
            
        end
        
        if break_out_of_loop
            break
        end
        
    end
    
    for nn = 1:sz(3)
        cur_data(:,:,nn,nnn) = mprfShiftTimeSeries(cur_data(:,:,nn,nnn),1000,'align_to_phase',cur_ph,10,false);
        
    end
    
end


set(0,'CurrentFigure',fh);

text_xpos = (size(cur_data,1) .* (0:size(cur_data,2)-1)) + size(cur_data,1)/2;



if size(cur_data,2) > 1
    cur_data = reshape(cur_data,[size(cur_data,1) .* size(cur_data,2), size(cur_data,3), size(cur_data,4)]);
else
    cur_data = squeeze(cur_data);
end

for n = 1:size(cur_data,3)
    subplot(size(cur_data,3),1,n);
    plot(cur_data(:,:,n));
    title(sprintf('Channel nr: %d',plot_data.channels(n)));
    ypos = max(max(cur_data(:,:,n)));
    text(text_xpos',ypos(:,ones(1,size(text_xpos,2)))',num2str(plot_data.periods'))
    
    if n == size(cur_data,3)
        xlabel('Samples')
        
    end
    
    ylabel('Tesla')
    
end

plot_data.cur_data = cur_data;
plot_data.channels = plot_data.channels;
plot_data.periods = plot_data.periods;
plot_data.repeats = plot_data.repeats;

set(fh, 'UserData',plot_data);





end

function plot_tseries(hobj, junk)

global mprf_fh
global data
params = get(mprf_fh,'UserData');

cur_periods = 1:size(data.data,2);


stim_freq = 10;
bl_freq = [9 11];
freq_range = sort([stim_freq bl_freq]);

nb_idx_tmp = round(mprfFreq2Index(size(data.data,1),bl_freq,1000));
sl_idx_tmp = round(mprfFreq2Index(size(data.data,1),stim_freq,1000));

snr_tmp = nan(length(cur_periods), size(params.plot.channels,2));
% snr_spec = nan(size(snr_tmp));
snr_data = nan([size(snr_tmp) size(data.data,3)]);

for n = 1:length(cur_periods)
    for nn = 1:size(params.plot.channels, 2);
        
        cur_data = squeeze(data.data(:,cur_periods(n),: , params.plot.channels(nn)));
        
        if any(~isnan(cur_data(:)))
            data_av_tmp = nanmean(cur_data,2);
            data_av_spec = cur_data(:);
            %             data_av_spec = data_av_spec(~isnan(data_av_spec));
            
            
            %             idx_range_spec = round(mprfFreq2Index(size(data_av_spec,1),freq_range,1000));
            %             nb_idx_spec = [idx_range_spec(1):idx_range_spec(2) idx_range_spec(3) : idx_range_spec(4)];
            %             sl_idx_spec = idx_range_spec(2) : idx_range_spec(3);
            
            ft_av_tmp = fft(data_av_tmp);
            %             ft_av_spec = fft(data_av_spec);
            ft_data = fft(cur_data);
            
            amp_av_tmp = (2 .* abs(ft_av_tmp) ./ size(ft_av_tmp,1));
            %             amp_av_spec = (2 .* abs(ft_av_spec) ./ size(ft_av_spec,1));
            amp_data = (2 .* abs(ft_data) ./ size(ft_data,1));
            
            snr_tmp(n,nn) = exp(nanmean(log(amp_av_tmp(sl_idx_tmp)))) ...
                ./ exp(nanmean(log(amp_av_tmp(nb_idx_tmp))));
            
            
            %             snr_spec(n,nn) = exp(nanmean(log(amp_av_spec(sl_idx_spec)))) ...
            %                 ./ exp(nanmean(log(amp_av_spec(nb_idx_spec))));
            
            snr_data(n,nn,:) = exp(nanmean(log(amp_data(sl_idx_tmp,:)))) ...
                ./ exp(nanmean(log(amp_data(nb_idx_tmp,:))));
            
        end
    end
end

av_amp = nanmean(snr_data,3);
ste_amp = nanstd(snr_data,[],3)./ sqrt(size(snr_data,3));

y_lims = [0.75 round(max(av_amp(:)))];
y_lims_02 = [0.75 round(max(snr_tmp(:)))];


for n = 1:size(params.plot.channels, 2);
    
    
    figure; hold on;
    
    p_stim_periods = [reshape(params.periods.stim,[],5); ...
        flipud(reshape(params.periods.stim,[],5))];
    
    stim_y = y_lims(ones(size(p_stim_periods,1)/2,1),:);
    stim_y = stim_y(:);
    
    p_blink_periods = [reshape(params.periods.blink,[],6); ...
        flipud(reshape(params.periods.blink,[],6))];
    
    
    blink_y = y_lims(ones(size(p_blink_periods,1)/2,1),:);
    blink_y = blink_y(:);
    
    p_blank_periods = [reshape(params.periods.blank,[],6); ...
        flipud(reshape(params.periods.blank,[],6))];
    
    blank_y = y_lims(ones(size(p_blank_periods,1)/2,1),:);
    blank_y = blank_y(:);
    
    
    patch(p_stim_periods,stim_y(:,ones(1,size(p_stim_periods,2))),'g')
    patch(p_blink_periods,blink_y(:,ones(1,size(p_blink_periods,2))),'r')
    patch(p_blank_periods,blank_y(:,ones(1,size(p_blank_periods,2))),'y')
    
    plot(cur_periods, av_amp(:,n),'k-','LineWidth',4);
    cur_lim = axis;
    errorbar(cur_periods,av_amp(:,n), ste_amp(:,n),'k.')
    plot(cur_periods,ones(1,size(av_amp,1)),'k--')
    axis(cur_lim);
    ylim(y_lims);
    title(sprintf('Channel nr: %d',params.plot.channels(n)));
    
    figure; hold on;
    
    p_stim_periods = [reshape(params.periods.stim,[],5); ...
        flipud(reshape(params.periods.stim,[],5))];
    
    stim_y = y_lims_02(ones(size(p_stim_periods,1)/2,1),:);
    stim_y = stim_y(:);
    
    p_blink_periods = [reshape(params.periods.blink,[],6); ...
        flipud(reshape(params.periods.blink,[],6))];
    
    
    blink_y = y_lims_02(ones(size(p_blink_periods,1)/2,1),:);
    blink_y = blink_y(:);
    
    p_blank_periods = [reshape(params.periods.blank,[],6); ...
        flipud(reshape(params.periods.blank,[],6))];
    
    blank_y = y_lims_02(ones(size(p_blank_periods,1)/2,1),:);
    blank_y = blank_y(:);
    
    
    patch(p_stim_periods,stim_y(:,ones(1,size(p_stim_periods,2))),'g')
    patch(p_blink_periods,blink_y(:,ones(1,size(p_blink_periods,2))),'r')
    patch(p_blank_periods,blank_y(:,ones(1,size(p_blank_periods,2))),'y')
    
    
    
    plot(cur_periods, snr_tmp(:,n),'k-','LineWidth',4);
    cur_lim = axis;
    plot(cur_periods,ones(1,size(snr_tmp,1)),'k--')
    axis(cur_lim);
    ylim(y_lims_02);
    title(sprintf('Channel nr: %d',params.plot.channels(n)));
    %
    %     figure; hold on;
    %
    %     p_stim_periods = [reshape(params.periods.stim,[],5); ...
    %         flipud(reshape(params.periods.stim,[],5))];
    %
    %     stim_y = y_lims(ones(size(p_stim_periods,1)/2,1),:);
    %     stim_y = stim_y(:);
    %
    %     p_blink_periods = [reshape(params.periods.blink,[],6); ...
    %         flipud(reshape(params.periods.blink,[],6))];
    %
    %
    %     blink_y = y_lims(ones(size(p_blink_periods,1)/2,1),:);
    %     blink_y = blink_y(:);
    %
    %     p_blank_periods = [reshape(params.periods.blank,[],6); ...
    %         flipud(reshape(params.periods.blank,[],6))];
    %
    %     blank_y = y_lims(ones(size(p_blank_periods,1)/2,1),:);
    %     blank_y = blank_y(:);
    %
    %
    %     patch(p_stim_periods,stim_y(:,ones(1,size(p_stim_periods,2))),'g')
    %     patch(p_blink_periods,blink_y(:,ones(1,size(p_blink_periods,2))),'r')
    %     patch(p_blank_periods,blank_y(:,ones(1,size(p_blank_periods,2))),'y')
    %
    %
    %     plot(cur_periods, snr_spec(:,n),'k-','LineWidth',4);
    %     cur_lim = axis;
    %     plot(cur_periods,ones(1,size(snr_spec,1)),'k--')
    %     axis(cur_lim);
    %     ylim(y_lims);
    %     title(sprintf('Channel nr: %d',params.plot.channels(n)));
    %
    
end
end
% function plot_amp_tseries(hobj, junk)
%
% global data
% global mprf_fh
%
% params = get(mprf_fh,'UserData');
% source_dir = fileparts(params.paths.source_file);
% cur_dir = pwd;
% cd(source_dir);
% make_gui_images = true;
% if isfield(params,'paths')
%     if isfield(params.paths,'amp_head_pp_file')
%
%         if isempty(params.paths.amp_head_pp_file)
%             [fname, fpath] = uigetfile('Please select file to plot');
%
%             if fname
%                 load(fullfile(fpath,fname));
%                 make_gui_images = false;
%
%
%             end
%         else
%
%             try
%                 load(params.paths.amp_head_pp_file)
%
%                 make_gui_images = false;
%
%             catch
%
%                 fprintf('Could not load stored images, creating new ones\n');
%
%             end
%         end
%
%
%     end
% end
%
% cd(cur_dir)
%
% if make_gui_images
%     stimfiles = dir(fullfile(params.paths.stim_dir,'*.mat'));
%     stim = load(fullfile(params.paths.stim_dir, stimfiles(1).name));
%
%     clear stimfiles
%
%     stim_range = double([min(stim.stimulus.images(:)) max(stim.stimulus.images(:))]);
%     bk = double(mode(mode(stim.stimulus.images(:,:,1))));
%     new_period = stim.stimulus.seq(stim.stimulus.trigSeq > 0);
%
%     im_out_size = [101 101];
%     im_seq = zeros([im_out_size, length(new_period)],'uint8');
%
%
%     for n = 1:length(new_period)
%         tmp_im = double(stim.stimulus.images(:,:,new_period(n)));
%         tmp_im = imresize(tmp_im,im_out_size,'nearest');
%         im_seq(:,:,n) = uint8(ceil(abs(tmp_im - bk) ./ max(abs(stim_range - bk))));
%
%     end
%
%     cur_data = data.data(:,params.periods.stim,:,:);
%     cur_data = reshape(cur_data,[],size(data.data,4));
%     fprintf('Computing amplitude modulation for %d channels:\n', size(cur_data,2))
%
%     stim_freq = 10;
%     bl_freq = [9 11];
%
%     snr = ones(size(data.data,2),size(data.data,4));
%     %tmp_data = squeeze(nanmean(data.data,3));
%
%     for n = 1:size(data.data,4)
%         fprintf('%d.',n);
%
%         tmp_data_02 = data.data(:,:,:,n);
%
%         if any(~isnan(tmp_data_02(:)))
%             tmp = fft(tmp_data_02,[],1);
%             tmp = tmp(1:floor(size(tmp,1)/2)+1,:,:);
%
%             nb_idx = round(mprfFreq2Index(size(tmp_data_02,1),bl_freq,1000));
%             sl_idx = round(mprfFreq2Index(size(tmp_data_02,1),stim_freq,1000));
%
%             spec_amp_nb = abs(tmp(nb_idx,:,:));
%             spec_amp_sl = abs(tmp(sl_idx,:,:));
%
%             av_nb = squeeze(exp(nanmean(log(spec_amp_nb))));
%             av_sl = squeeze(spec_amp_sl);
%
%             snr(:,n) = nanmedian(av_sl,2) ./ nanmedian(av_nb,2);
%             %    fprintf('Amplitude modulation for channel %d: %0.2f\n',...
%             %        plot_data.channels(nn), snr)
%
%         end
%
%     end
%
%
%
%     fprintf('Done\n')
%
%
%     fh_02 = figure;
%
%
%
%     for n = 1:size(snr,1)
%
%         fh_02 = megPlotMap(snr(n,:),[0.75 max(snr(:)./0.05) .* 0.05],fh_02,jet(256),sprintf('Stimulus position %d',n));
%         drawnow;
%
%         tmp  = getframe(fh_02);
%
%         if n == 1
%             head_im = zeros([size(tmp.cdata) size(snr,1)],'uint8');
%
%         end
%         head_im(:,:,:,n) = tmp.cdata;
%
%
%     end
%
%     cur_time = datestr(now);
%     cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
%
%     source_dir = fileparts(params.paths.source_file);
%
%     save(fullfile(source_dir,['Amplitude_head_overlay_per_period_' cur_time]),'head_im','im_seq')
%
%
%     params.paths.amp_head_pp_file = fullfile(params.paths.data_dir,['Amplitude_head_overlay_per_period_' cur_time]);
%
% end
%
%
%
%
% mprfSNRPerPeriodGUI(params, im_seq,head_im);
%
%
% end


%
%
% function make_phase_figure(hobj, junk)
%
% global mprf_fh
% global data
%
% fh = figure; hold on;
% if isobject(fh)
%     cur_fh = str2double(inputdlg('Which figure?','',1,{num2str(fh.Number-1)}));
%
% else
%     cur_fh = str2double(inputdlg('Which figure?','',1,{num2str(fh-1)}));
%
% end
%
%
% plot_data = get(cur_fh,'UserData');
% % params = get(mprf_fh,'UserData');
%
% cur_data = reshape(plot_data.cur_data,[],numel(plot_data.periods), numel(plot_data.repeats), numel(plot_data.channels));
%
%
% if size(cur_data,3) > 1
%     cur_data_02 = reshape(cur_data,size(cur_data,1), [], size(cur_data,4));
% end
%
%
% n_plots = size(cur_data_02,3);
%
% while numel(factor(n_plots)) == 1;
%
%     n_plots = n_plots + 1;
% end
%
%
% tmp = factor(n_plots);
% nrows = prod(tmp(2:2:length(tmp)));
% ncols = prod(tmp(1:2:length(tmp)));
%
%
% ft_data = fft(cur_data_02,[],1);
%
% ph = angle(ft_data(mprfFreq2Index(size(cur_data,1),10,1000),:,:));
%
%
% for n = 1:size(ph,3)
%     subplot(nrows,ncols,n);
%     [Y, X] = hist(ph(1,:,n),-pi:pi/5:pi);
%     bar(X,Y,'FaceColor',[0.5 0.5 0.5],...
%         'EdgeColor','k')
%     hold on;
%     tmp_y = ylim;
%     tmp_x = prctile(ph(1,:,n),[50 25 75]);
%     plot(tmp_x([1 1]), tmp_y,'k-','LineWidth',2)
%     plot(tmp_x([2 3; 2 3]), ylim, 'k--','LineWidth',2)
%     title(sprintf('Channel %d',plot_data.channels(n)));
%     hold off
%
%
%
% end
%
%
%
%
% end


%
% function make_amp_head_layout(hobj, junk)
%
%
%
% global data
% global mprf_fh
%
% params = get(mprf_fh,'UserData');
%
% bl_range = 5;
% sl_range = 1;
% stim_freq = 10;
% bl_start_end = stim_freq + [-bl_range/2 bl_range/2];
% sl_start_end = stim_freq + [-sl_range/2 sl_range/2];
%
% freq_range = [bl_start_end(1) sl_start_end bl_start_end(2)];
%
% fprintf('Baseline range: %0.2f : %0.2f\n',bl_start_end)
% fprintf('Stimulus range: %0.2f : %0.2f\n',sl_start_end)
%
% cur_data = data.data(:,params.periods.stim,:,:);
% cur_data = reshape(cur_data,[],size(data.data,4));
% snr = ones(1,size(cur_data,2));
% fprintf('Computing amplitude modulation for %d channels:\n', size(cur_data,2))
% for n = 1:size(cur_data,2)
%     fprintf('%d.',n);
%     tmp_data = cur_data(~isnan(cur_data(:,n)),n);
%
%     if any(tmp_data)
%         tmp = fft(tmp_data);
%         tmp = tmp(1:floor(size(tmp,1)/2)+1,:,:);
%
%
%         idx_range = round(mprfFreq2Index(size(tmp_data,1),freq_range,1000));
%
%         nb_idx = [idx_range(1):idx_range(2) idx_range(3) : idx_range(4)];
%         sl_idx = idx_range(2) : idx_range(3);
%
%         spec_amp_nb = abs(tmp(nb_idx));
%         spec_amp_sl = abs(tmp(sl_idx));
%
%         av_nb = exp(nanmean(log(spec_amp_nb)));
%         av_sl = exp(nanmean(log(spec_amp_sl)));
%
%         snr(n) = av_sl ./ av_nb;
%         %    fprintf('Amplitude modulation for channel %d: %0.2f\n',...
%         %        plot_data.channels(nn), snr)
%
%     end
%
% end
% fprintf('Done\n')
% fprintf('Making head layout plots...\n')
%
% fh = figure;
% megPlotMap(snr,[1 ceil(max(snr ./ 0.05)) .* 0.05],fh,jet(256),'SNR');
%
%
% get_layout_data = true;
% if isfield(params,'plot')
%     if isfield(params.plot,'head')
%         if isfield(params.plot.head,'xpos')
%             get_layout_data = false;
%
%         end
%     end
% end
%
% if get_layout_data
%     load(params.paths.example_hdr)
%     layout = ft_prepare_layout(hdr);
%
%     params.plot.head.xpos = layout.pos(1:157,1);
%     params.plot.head.ypos = layout.pos(1:157,2);
%     params.plot.head.outline = layout.outline;
% end
%
%
% cur_chan = 1:157;
%
%
% figure;
% hold on
% for n = 1:length(params.plot.head.outline)
%     plot(params.plot.head.outline{n}(:,1), params.plot.head.outline{n}(:,2), 'k-','LineWidth',2)
% end
%
% plot(params.plot.head.xpos(1:157),params.plot.head.ypos(1:157),'k.')
% text(params.plot.head.xpos(cur_chan)',params.plot.head.ypos(cur_chan)',num2str(cur_chan'))
% axis square
% axis off
%
% set(mprf_fh,'UserData',params);
%
%
% end

