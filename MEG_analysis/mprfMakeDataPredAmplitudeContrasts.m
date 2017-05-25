function  mprfMakeDataPredAmplitudeContrasts(data, pred, stim, periods)

im_seq = uint8(stim);
%
% params = get(mprf_fh,'UserData');
%
% stimfiles = dir(fullfile(params.paths.stim_dir,'*.mat'));
% stim = load(fullfile(params.paths.stim_dir, stimfiles(1).name));
%
% clear stimfiles
%
% stim_range = double([min(stim.stimulus.images(:)) max(stim.stimulus.images(:))]);
% bk = double(mode(mode(stim.stimulus.images(:,:,1))));
% new_period = stim.stimulus.seq(stim.stimulus.trigSeq > 0);
%
% im_out_size = [101 101];
% im_seq = zeros([im_out_size, length(new_period)],'uint8');
%
%
% for n = 1:length(new_period)
%     tmp_im = double(stim.stimulus.images(:,:,new_period(n)));
%     tmp_im = imresize(tmp_im,im_out_size,'nearest');
%     im_seq(:,:,n) = uint8(ceil(abs(tmp_im - bk) ./ max(abs(stim_range - bk))));
%
% end

blank_periods = periods.blank;
blink_periods = periods.blink;
stim_periods = periods.stim;

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

snr_temp_data = ones(size(data.data,2),size(data.data,4));
snr_temp_pred = ones(size(data.data,2),size(data.data,4));

for n = 1:size(data.data,4)
    fprintf('%d.',n);
    
    for nn = 1:2
        if nn == 1
            tmp_data = data.data(:,:,:,n);
            
        elseif nn == 2
            tmp_data = pred.data(:,:,:,n);
            
        end
        
        if any(~isnan(tmp_data(:)))
            tmp = fft(tmp_data,[],1);
            tmp = tmp(1:floor(size(tmp,1)/2)+1,:,:);
            
            nb_idx = round(mprfFreq2Index(size(tmp_data,1),bl_freq,1000));
            sl_idx = round(mprfFreq2Index(size(tmp_data,1),stim_freq,1000));
            
            spec_amp_nb = abs(tmp(nb_idx,:,:));
            spec_amp_sl = abs(tmp(sl_idx,:,:));
            
            av_nb = squeeze(exp(nanmean(log(spec_amp_nb))));
            av_sl = squeeze(spec_amp_sl);
            
            if nn == 1
                snr_temp_data(:,n) = nanmedian(av_sl,2) ./ nanmedian(av_nb,2);
            elseif nn == 2
                snr_temp_pred(:,n) = nanmedian(av_sl,2) ./ nanmedian(av_nb,2);
                
            end
            
            
            
        end
        
    end
    
    
    
end



fprintf('Done\n')
snr_temp_data(isnan(snr_temp_data)) = 0;
snr_temp_pred(isnan(snr_temp_pred)) = 0;

for nn = 1:2
    if nn == 1
        snr_temp = snr_temp_data;
        title_str = 'Measured';
        ui_pos = 'northwest';
        
    elseif nn == 2
        snr_temp = snr_temp_pred;
        title_str = 'Predicted';
        ui_pos = 'northeast';
    end
    
    
    stim_vs_blank_temp = stim_vs_blank * snr_temp;
    lhs_vs_rhs_temp = lhs_vs_rhs * snr_temp;
    up_vs_low_temp = up_vs_low * snr_temp;
    
    up_vs_blank_temp = up_vs_blank * snr_temp;
    low_vs_blank_temp = low_vs_blank * snr_temp;
    left_vs_blank_temp = left_vs_blank * snr_temp;
    right_vs_blank_temp = right_vs_blank * snr_temp;
    
    
    f = figure;
    megPlotMap(stim_vs_blank_temp,[floor(min(stim_vs_blank_temp./0.05) .* 0.05) ceil(max(stim_vs_blank_temp./0.05) .* 0.05)],...
        f,jet(256), ['Stim vs blank,' title_str]);
    movegui(f, ui_pos);
    
    f = figure;
    megPlotMap(lhs_vs_rhs_temp,[floor(min(left_vs_blank_temp./0.05) .* 0.05) ceil(max(left_vs_blank_temp./0.05) .* 0.05)], ...
        f,jet(256), ['left vs right,' title_str]);
    movegui(f, ui_pos);
    
    f = figure;
    megPlotMap(up_vs_low_temp,[floor(min(up_vs_low_temp./0.05) .* 0.05) ceil(max(up_vs_low_temp./0.05) .* 0.05)], ...
        f,jet(256), ['up vs low,' title_str]);
    movegui(f, ui_pos);
    
    
    f = figure;
    megPlotMap(up_vs_blank_temp,[floor(min(up_vs_blank_temp./0.05) .* 0.05) ceil(max(up_vs_blank_temp./0.05) .* 0.05)],...
        f,jet(256), ['Up vs blank,' title_str]);
    movegui(f, ui_pos);
    
    
    f = figure;
    megPlotMap(low_vs_blank_temp,[floor(min(low_vs_blank_temp./0.05) .* 0.05) ceil(max(low_vs_blank_temp./0.05) .* 0.05)],...
        f,jet(256), ['Low vs blank,' title_str]);
    movegui(f, ui_pos);
    
    
    f = figure;
    megPlotMap(left_vs_blank_temp,[floor(min(left_vs_blank_temp./0.05) .* 0.05) ceil(max(left_vs_blank_temp./0.05) .* 0.05)], ...
        f,jet(256), ['left vs blank,' title_str]);
    movegui(f, ui_pos);
    
    f = figure;
    megPlotMap(right_vs_blank_temp,[floor(min(right_vs_blank_temp./0.05) .* 0.05) ceil(max(right_vs_blank_temp./0.05) .* 0.05)], ...
        f,jet(256), ['right vs blank,' title_str]);
    movegui(f, ui_pos);
    
    
    
end



