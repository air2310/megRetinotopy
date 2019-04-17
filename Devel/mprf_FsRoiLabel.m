
fs_roi_vz = '/mnt/storage_2/MEG/Retinotopy/Quality_check/wlsubj004/rois/surface/freesurfer/FS_visualize';

surfaces_to_load = {'lh.white','rh.white'};

for n = 1:length(surfaces_to_load)
    
    cur_surf = fullfile(FS_surface,surfaces_to_load{n});
    
    tmp = strsplit(surfaces_to_load{n},'.');
    cur_hs = tmp{1};
    fprintf('Exporting parameters for %s hemisphere:\n',cur_hs);
    
    
    % load our rois which are now as morphometry file formats (saved using
    % write_curv function)
    roi_name = roi_names;
    
    for roi_idx = 1:length(roi_names)
        
        cur_roi = roi_name{roi_idx};
        
        [~,r_name] = fileparts(roi_names{n});
        
        
        
        tmp_roi = read_curv(fullfile(fs_roi_dir,strcat(cur_hs,'.',cur_roi(2:end-4))));
        lindex = find(~isnan(tmp_roi));
        
        %create a label file on the terminal as lh.V1.label for eg:
        labelfile = fullfile(fs_roi_vz,strcat(cur_hs,'.',cur_roi(2:end-4),'.label')); % 'lh.V1.label';
        
        npoints1 = length(lindex);
        npoints = npoints1;
        lxyz = zeros(npoints,3);
        lvals  = zeros(npoints,1);
        
        % open as an ascii file
        fid = fopen(labelfile, 'w') ;
        if(fid == -1)
            fprintf('ERROR: could not open %s\n',labelfile);
            return;
        end
        
        
        % Make sure they are npoints by 1 %
        lindex = reshape(lindex,[npoints 1]);
        lxyz   = reshape(lxyz,[npoints 3]);
        lvals  = reshape(lvals,[npoints 1]);
        
        l = [lindex lxyz lvals];
        fprintf(fid,'%d %f %f %f %f\n',l') ;
        
        fclose(fid) ;
    end
end


%% Saving prf data for visualizing on freesurfer surface
% prf_data folder containing prf parameters in different formats. We need
% the prf parameters in the freesurfer surface in order to run this section
% of the code.

sub = 'wlsubj068';
FS_surface = sprintf('/mnt/storage_2/MEG/Retinotopy/Data/Freesurfer_directory/%s/surf',sub);
fs_prf_data = sprintf('/mnt/storage_2/MEG/Retinotopy/Subject_sessions/%s/prf_data/surface/freesurfer',sub);
fs_prf_data_vz = sprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/prf_data/surface/freesurfer/FS_visualize',sub);
   
if ~exist(fs_prf_data_vz,'dir')
    mkdir(fs_prf_data_vz)
end

surfaces_to_load = {'lh.white','rh.white'};

for n = 1:length(surfaces_to_load)
    
    cur_surf = fullfile(FS_surface,surfaces_to_load{n});
    
    tmp = strsplit(surfaces_to_load{n},'.');
    cur_hs = tmp{1};
    fprintf('Exporting parameters for %s hemisphere:\n',cur_hs);
   
    
    %lh_files = dir(fullfile(fs_prf_data,'lh.*'));
    h_files = dir(strcat(fs_prf_data,'/',cur_hs,'.*'));
    
    for nn = 1:length(h_files)
        
        cur_h_file = h_files(nn).name;
        [~,~,par_name] = fileparts(cur_h_file);
        tmp_pm = cur_h_file; %'lh.polar_angle';
        
        tmp_pm_fpath = fullfile(fs_prf_data,tmp_pm);
        tmp_pm_file = read_curv(tmp_pm_fpath);
        idx_nan = find(isnan(tmp_pm_file));
        tmp_pm_file(idx_nan) = 0;
        
        write_curv(fullfile(fs_prf_data_vz,tmp_pm),tmp_pm_file,0);
        
        fprintf('\n Max of %s = %f \n ',cur_h_file,max(tmp_pm_file));
        fprintf('\n Min of %s = %f \n ',cur_h_file,min(tmp_pm_file));
        
    end
end

% Once this is run, visualize the prf parameters on freesurfer surface
% using freeview
% We need to have lh.inflated and rh.inflated freesurfer surface
% freeview -f lh.inflated (in the terminal)
% Load the prf data (eccentricity smoothed, polar angle smoothed and betas smoothed for now) as overlay on the surface



%% CHECKING SOME STUFF (delete later)
read_lbl = 0;
if read_lbl == 1
        
    fpath = '/mnt/storage_2/MEG/Retinotopy/Quality_check/wlsubj004/rois/surface/freesurfer/FS_visualize';
    fname = 'lh.V1_check.label';
    roi_file = fullfile(fpath,fname);
    
    % open it as an ascii file
    fid = fopen(roi_file, 'r') ;
    if(fid == -1)
        fprintf('ERROR: could not open %s\n',fname);
        return;
    end
    
    fgets(fid);
    if(fid == -1)
        fprintf('ERROR: could not open %s\n',fname);
        return;
    end
    
    line = fgets(fid) ;
    nv = sscanf(line, '%d') ;
    l = fscanf(fid, '%d %f %f %f %f\n') ;
    l = reshape(l, 5, nv) ;
    l = l' ;
    
    fclose(fid) ;    
   
    
    labelfile = '/mnt/storage_2/MEG/Retinotopy/Quality_check/wlsubj004/rois/surface/freesurfer/FS_visualize/lh.V1_check.label';
      % open as an ascii file
        fid = fopen(labelfile, 'w') ;
        if(fid == -1)
            fprintf('ERROR: could not open %s\n',labelfile);
            return;
        end
        
        lindex = l(:,1);
        npoints1 = length(lindex);
        npoints = npoints1;
        lxyz = zeros(npoints,3);
        lvals  = zeros(npoints,1);
        
        % Make sure they are npoints by 1 %
        lindex = reshape(lindex,[npoints 1]);
        lxyz   = reshape(lxyz,[npoints 3]);
        lvals  = reshape(lvals,[npoints 1]);
        
        l = [lindex lxyz lvals];
        fprintf(fid,'%d %f %f %f %f\n',l') ;
        
        fclose(fid) ;
    
    
    
    
    
end


%% 19 reference phases for two channels with high VE and two with low VE
clf;
sub_name = 'wlsubj040';

main_dir = sprintf('/mnt/storage_2/MEG/Retinotopy/Subject_sessions/%s/modeling/results/original_model/',sub_name);
save_dir = sprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/Ph_opt',sub_name);
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end    

switch sub_name
    case 'wlsubj004'
        Results_path = strcat(main_dir,'Run_Stimulus_locked_model_fit_lo_15_Feb_2019_17_25_02');
    case 'wlsubj030'
        Results_path = strcat(main_dir,'Run_Stimulus_locked_model_fit_lo_13_Feb_2019_12_06_47');
    case 'wlsubj040'
        Results_path = strcat(main_dir,'Run_Stimulus_locked_model_fit_lo_13_Feb_2019_10_09_56');
    case 'wlsubj058'
        Results_path = strcat(main_dir,'Run_Stimulus_locked_model_fit_lo_12_Feb_2019_16_33_44');
    case 'wlsubj068'
        Results_path = strcat(main_dir,'Run_Stimulus_locked_model_fit_lo_12_Feb_2019_16_58_34');
end
        
% load results from mprfSession_run_model_server.m
%Results_path = '/mnt/storage_2/MEG/Retinotopy/Subject_sessions/wlsubj058/modeling/results/original_model/Run_Stimulus_locked_model_fit_lo_12_Feb_2019_16_33_44';
results_f = fullfile(Results_path,'Results.mat');
load(results_f);

% Rank channels in increasing order of mean variance explained from the 19
% leave one out iterations.

% Use this to decide the order in which images are plotted in the subplots
% Low ve on the top - mostly in the frontal areas
% high VE bottom - occipital areas
PH_opt = results.PH_opt;

[ve_sort, chan_sort] = sort(results.corr_mat);

chan = [chan_sort(1:3) chan_sort(sum(~isnan(ve_sort))/2:sum(~isnan(ve_sort))/2+1)  chan_sort(sum(~isnan(ve_sort))-2:sum(~isnan(ve_sort)))] ;


% idx_high = find(mean(results.VE_opt,1) == max(mean(results.VE_opt,1)));
% idx_low = find(mean(results.VE_opt,1) == min(mean(results.VE_opt,1)));
%idx_high = 1;


% Highlight the channels on the head plot
figPoint_rlchve_map = figure(100);
subplot(3,3,5);
hold_figure.flag =1;
hold_figure.figh = figPoint_rlchve_map;
figPoint_rlchve_map = mprfPlotHeadLayout(chan,0,5,0,hold_figure);

% interpmethod = [];
% subplot(3,3,5);
% megPlotMap(results.corr_mat,[0 0.6],figPoint_rlchve_map,'jet',...
%     'Phase ref fit stimulus locked',[],[],'interpmethod',interpmethod);

plot_loc = setdiff(1:length(chan),5);
plot_loc = [plot_loc plot_loc(end)+1];

for idx_chan = 1:length(chan)
    
cur_chan = chan(idx_chan);

idx_plot = plot_loc(idx_chan);

subplot(3,3,idx_plot);    

%figure;

mprf_polarplot(ones(size(PH_opt,1),1),PH_opt(:,cur_chan));
x_bl = ones(size(PH_opt(:,cur_chan))).* cos(PH_opt(:,cur_chan));
y_bl = ones(size(PH_opt(:,cur_chan))).* sin(PH_opt(:,cur_chan));
for ct = 1:size(PH_opt(:,cur_chan),1)
    hold on; plot([0 x_bl(ct)],[0 y_bl(ct)],'ro-');
end

chan_title = sprintf('CH: %d \n',cur_chan );
ve_title = sprintf('VE: %f', results.corr_mat(:,cur_chan)); 
%ph_title = sprintf('PH: %f \t', unique(PH_opt(:,cur_chan))');

%tle = sprintf(strcat(ve_title,'\n',ph_title));
tle = sprintf(strcat(chan_title,'\n',ve_title));
title(tle);

end

cd(save_dir)
saveas(figPoint_rlchve_map,strcat('ref_ph_19lo','.tif')); % original map


%% Plot the MEG scalp maps and range curves for every scaling values for pRF size and pRF position scaling
clear all;
close all;
sub_name = 'wlsubj004';

main_dir = sprintf('/mnt/storage_2/MEG/Retinotopy/Subject_sessions/%s/modeling/results/prf_size_range/',sub_name);
save_dir = sprintf('/mnt/storage_2/MEG/Retinotopy/Quality_check/%s/prf_size_range',sub_name);
if ~exist(save_dir,'dir')
    mkdir(save_dir);
end    

switch sub_name
    case 'wlsubj004'
        Results_path = strcat(main_dir,'Run_Stimulus_locked_model_fit_lo_15_Feb_2019_17_24_09');
    case 'wlsubj030'
        Results_path = strcat(main_dir,'Run_Stimulus_locked_model_fit_lo_13_Feb_2019_12_07_29');
    case 'wlsubj040'
        Results_path = strcat(main_dir,'Run_Stimulus_locked_model_fit_lo_13_Feb_2019_10_08_01');
    case 'wlsubj058'
        Results_path = strcat(main_dir,'Run_Stimulus_locked_model_fit_lo_12_Feb_2019_16_53_30');
    case 'wlsubj068'
        Results_path = strcat(main_dir,'Run_Stimulus_locked_model_fit_lo_12_Feb_2019_16_59_20');
end
        
% load results from mprfSession_run_model_server.m
%Results_path = '/mnt/storage_2/MEG/Retinotopy/Subject_sessions/wlsubj058/modeling/results/original_model/Run_Stimulus_locked_model_fit_lo_12_Feb_2019_16_33_44';
results_f = fullfile(Results_path,'Results.mat');
load(results_f);

% Size or position scaling iteration values
if strcmpi(model.type,'pRF size range')
    par_it = model.params.sigma_range;
elseif strcmpi(model.type,'position (x,y) range')
    par_it = model.params.x0_range;
end

% variance explained for each range iteration and channel as a 3-D matrix-
% 1 x # range iteration x # channels
corr_tmp = results.corr_mat;

% number of range iterations
n_par_it = size(corr_tmp,2);

% Selection of channels for averaging 
chan_sel = 'all';
if strcmpi(chan_sel,'back')
    load(which('meg160_example_hdr.mat'))
    layout = ft_prepare_layout([],hdr);
    xpos = layout.pos(1:157,1);
    ypos = layout.pos(1:157,2);
    %good_chan = find(ypos<1 & xpos<1);
    good_chan = find(ypos<0 & xpos<1);
    %n_iter = size(results.scrambled_corr,2);
    
elseif strcmpi(chan_sel,'rel')
    [rel_fname, rel_fpath] = uigetfile('*.mat','Select reliability run results');
    rel_results = load(fullfile(rel_fpath, rel_fname));
    med_corr = (nanmedian(rel_results.results.split_half.corr_mat,2));
    good_chan = find(med_corr>0.2);
    figPoint_occ_map = mprfPlotHeadLayout(good_chan);
    
elseif strcmpi(chan_sel,'all')
    load(which('meg160_example_hdr.mat'))
    layout = ft_prepare_layout([],hdr);
    xpos = layout.pos(1:157,1);
    ypos = layout.pos(1:157,2);
    %good_chan = find(ypos<1 & xpos<1);
    good_chan = find(ypos<0 | ypos>=0 | xpos<0 | xpos>=0);

end



% Figure for 
% graph for variance explained vs pRF size
chan_ave = 'False';
switch chan_ave
    case 'True'
        % Select only the channels from good_chan
        tmp_corr_sl = nan(n_par_it,length(good_chan));
        all_corr = results.corr_mat;
        for i=1:n_par_it
            tmp_corr_sl(i,:) = squeeze(all_corr(:,i,good_chan,:,1));
        end
        
        % Average variance explained values for all the selected channels
        corr_avg_sl = nanmean(tmp_corr_sl,2);
        corr_std_sl = nanstd(tmp_corr_sl,0,2);
        corr_ste_sl = corr_std_sl ./ sqrt(size(tmp_corr_sl,2));
        corr_CI_sl = 1.96 .* corr_ste_sl;
        
        fh_sl_VEparams = figure; set(gcf, 'Position', [1000, 592, 838, 746]);
        
        % Create shaded error bar with 'patch' function
        lo = corr_avg_sl - corr_CI_sl;
        hi = corr_avg_sl + corr_CI_sl;
        color = [0.5 0.5 0.5];
        err= patch([par_it, fliplr(par_it)], [lo', fliplr(hi')], color, 'FaceAlpha', 0.5, 'LineStyle',':');
        
        hold on;
        plot(par_it,corr_avg_sl,'r','Linewidth',3);
        
        % Set the axis
        set(gca,'TickDir', 'out');
        if strcmpi(model.type,'prf size range')
            xlabel('Size scaling factor (ratio)');
            set(gca,'XTick', [0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10],'XTickLabel',{'0.2', '0.5', '1.0', '1.5', '3.0', '5.0', '10'});
            set(gca, 'XScale', 'log');
        elseif strcmpi(model.type,'position (x,y) range')
            xlabel('position (deg)');
            set(gca,'XTick', par_it,'XTickLabel',rad2deg(model.params.x0_range));
        end
        xlim([par_it(1),par_it(end)]);
        title(sprintf('Variance explained %s sl locked',model.type));
        set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
        ylabel('Variance explained (%)');
        
    case 'False'
        
        % Select only the channels from good_chan
        tmp_corr_sl = nan(n_par_it,length(good_chan));
        all_corr = results.corr_mat;
        for i=1:n_par_it
            tmp_corr_sl(i,:) = squeeze(all_corr(:,i,good_chan,:,1));
        end
        
        % Average variance explained values for all the selected channels
        corr_sl = tmp_corr_sl;
        %corr_std_sl = nanstd(tmp_corr_sl,0,2);
        %corr_ste_sl = corr_std_sl ./ sqrt(size(tmp_corr_sl,2));
        %corr_CI_sl = 1.96 .* corr_ste_sl;
        
        peak_par_it = max(corr_sl,[],1);
        %figure, hist(peak_par_it);
        
        idx_peak_par_it_tmp = (corr_sl == repmat(peak_par_it,[size(corr_sl,1),1]));
        
        idx_peak_par_it_l = nan(1,size(corr_sl,2));
        idx_peak_par_it_h = nan(1,size(corr_sl,2));
        for idx_chan = 1:size(corr_sl,2)
            % There are some Nan values in the variance explained matrix,
            % which causes the peak values to be an exmpty matrix
            if sum(idx_peak_par_it_tmp(:,idx_chan))==1
                idx_peak_par_it_l(idx_chan) = ( par_it(7) <= par_it(find(idx_peak_par_it_tmp(:,idx_chan))) & par_it(find(idx_peak_par_it_tmp(:,idx_chan))) <= par_it(11));
                idx_peak_par_it_h(idx_chan) = ( par_it(7) > par_it(find(idx_peak_par_it_tmp(:,idx_chan))) | par_it(find(idx_peak_par_it_tmp(:,idx_chan))) > par_it(11));
            end
        end

        
        chan_peak_par_it_l = good_chan(find(idx_peak_par_it_l))';
        chan_peak_par_it_h = good_chan(find(idx_peak_par_it_h))';
        
        
        fh_sl_VEparams = figure; %set(gcf, 'Position', [1000, 592, 838, 746]);
        
        % Create shaded error bar with 'patch' function
        %lo = corr_avg_sl - corr_CI_sl;
        %hi = corr_avg_sl + corr_CI_sl;
        %color = [0.5 0.5 0.5];
        %err= patch([par_it, fliplr(par_it)], [lo', fliplr(hi')], color, 'FaceAlpha', 0.5, 'LineStyle',':');
        
        subplot(2,2,1);
        
        plot(par_it,squeeze(all_corr(:,:,chan_peak_par_it_h)),'g','Linewidth',0.5);
        hold on;
        plot(par_it,squeeze(all_corr(:,:,chan_peak_par_it_l)),'r','Linewidth',0.5);

        
        % Set the axis
        set(gca,'TickDir', 'out');
        if strcmpi(model.type,'prf size range')
            xlabel('Size scaling factor (ratio)');
            set(gca,'XTick', [0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10],'XTickLabel',{'0.2', '0.5', '1.0', '1.5', '3.0', '5.0', '10'});
            set(gca, 'XScale', 'log');
        elseif strcmpi(model.type,'position (x,y) range')
            xlabel('position (deg)');
            set(gca,'XTick', par_it,'XTickLabel',rad2deg(model.params.x0_range));
        end
        xlim([par_it(1),par_it(end)]);
        title(sprintf('Variance explained %s sl locked', model.type));
        set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
        ylabel('Variance explained (%)');
        
        
        % Make a subplot 2x2 and in the second section, plot the scalp map
        % with the electrodes in the corresponding colors
        %figPoint_rlchve_map = figure;
        subplot(2,2,2);
        hold_figure.flag =1;
        hold_figure.figh = fh_sl_VEparams;
        hold_figure.c = 'r';
        mprfPlotHeadLayout(chan_peak_par_it_l,0,10,0,hold_figure);
        
%       figPoint_rlchve_map = figure;
        subplot(2,2,3);
        hold_figure.flag =1;
        hold_figure.figh = fh_sl_VEparams;
        hold_figure.c = 'g';
        mprfPlotHeadLayout(chan_peak_par_it_h,0,10,0,hold_figure);
        
        %       figPoint_rlchve_map = figure;
        subplot(2,2,4);
        megPlotMap(results.corr_at_one_sl,[0 0.6],fh_sl_VEparams,'jet',...
            'Phase ref fit stimulus locked',[],[],'interpmethod',[]);
end

cd(save_dir)
range_type = model.type;
range_type(range_type == ' ') = '_';
saveas(fh_sl_VEparams,strcat(range_type,'all_Channels','.tif')); % original map




