function mprf_plot_figures(plot_type,interpmethod,compare_phfits)
% mprf_plot_figures - To make the figures from the results saved using
% mprfSession_run_model command
%
% Inputs -
%         plot_type : 1,2,3,4,5,6,7,8 for different figures as given below
%         interpmethod : 'nearest' for nearest interpolation
%
%  1 -  'Original model'
%  2 -  'scrambled'; channels for averaging  - all the channels on the back on head
%  3 -  'scrambled'; channels for averaging based on the reliability metric;
%  4 -  'scrambled individual channels' - p value for individual channels, to select the channels that are sensitive to scrambling;
%  5 -  'range'; channels for averaging  - all the channels in the back of head
%  6 -  'range'; channels for averaging based on the reliability metric;
%  7 -  'range'; channels for averaging are selected based on whether they are senstive to scrambling of pRF parameters;
%  8 -  'reliability'; - across repeat reliability metric
%  9 -  'phase_diff'; - Reliable phase value for the phase referenced
%                        amplitude fits

%*******************************************************************************************************************
% For original,  model_type = 'Original model'; chan_sel = '';
% For Scrambled, model_type = 'scrambled';
%                              chan_sel = 'rel' (for reliable channels);
%                              chan_sel = 'back'(for all electrodes in the back)
% For Range,     model_type = 'range';
%                              chan_sel = 'rel';
%                              chan_sel = 'back';
%                              chan_sel = 'ind_scr' (for channels with p value <0.01 for the scrambling condition)
% For reliability,
%                model_type = 'reliability';
% For phase_diff,model_type = 'phase_diff';
%*******************************************************************************************************************

%....................................................
% Original  -  Data_fit
%              model_fit
%
% Scrambled - Data_fit -  all scrambled -       back
%                                              rel
%                         group scrambled -     back
%                                              rel
%           - Model_fit - all scrambled        back
%                                              rel
%                         group scrambled      back
%                                              rel
% Range -     Data_fit - sigma                 back
%                                              rel
%                                              ind_scr
%                      - Position              back
%                                              rel
%                                              ind_scr
% Range -     Model_fit - sigma                 back
%                                              rel
%                                              ind_scr
%                      - Position              back
%                                              rel
%                                              ind_scr
%......................................................
global mprfSESSION
load mprfSESSION.mat;
try
    % if ~exist('plot_type','var') || isempty(plot_type)
    %     fprintf('No plot type selected, quitting\n');
    %     return;
    % end
    
    %if nargin > 1
    
    
    if ~exist('interpmethod','var') || isempty(interpmethod)
        interpmethod = [];
    end
    
    if ~exist('compare_phfits','var') || isempty(compare_phfits)
        compare_phfits = 0;
    end
    
    
    cur_time = mprf__get_cur_time;
    res_dir = mprf__get_directory('model_plots');
    %     main_dir = mprf__get_directory('main_dir');
    main_dir  = mprfSESSION.init.main_dir;
    
    %plot_type = 1;
    switch plot_type
        case 1
            model_type = 'Original model';  chan_sel = '';
            
        case 2
            model_type = 'scrambled'; chan_sel = 'back';
            
        case 3
            model_type = 'scrambled'; chan_sel = 'rel';
            
        case 4
            model_type = 'scrambled individual channel'; chan_sel = '';
            
        case 5
            model_type = 'range'; chan_sel = 'back';
            
        case 6
            model_type = 'range'; chan_sel = 'rel';
            
        case 7
            model_type = 'range'; chan_sel = 'ind_scr';
            
        case 8
            model_type = 'reliability'; chan_sel = '';
            
        case 9
            model_type = 'phase_diff'; chan_sel = '';
            
    end
    
    
    switch model_type
        case 'Original model'
            
            %% for Original
            % Select Original results
            [orig_fname, orig_fpath] = uigetfile('*.mat','Select original results');
            orig_results = load(fullfile(orig_fpath, orig_fname));
            
            % Original correlation
            fh_sl_map=figure;
            megPlotMap(orig_results.results.corr_mat,[0 0.6],fh_sl_map,'jet',...
                'Phase ref fit stimulus locked',[],[],'interpmethod',interpmethod);
            %saveas(fh_sl_map,strcat('sl_map','.jpg'));
            
            fh_sl_ph_map=figure;
            megPlotMap(rad2deg(orig_results.results.PH_opt(1,:)),[0 360],fh_sl_ph_map,'jet',...
                'Phase values (data driven)',[],[],'interpmethod',interpmethod);
            %saveas(fh_sl_ph_map,strcat('ph_map','.jpg'));
            
            chan = 15;
            figure,
            mprf_polarplot(ones(size(orig_results.results.PH_opt,1),1),orig_results.results.PH_opt(:,chan));
            x_bl = ones(size(orig_results.results.PH_opt(:,chan))).* cos(orig_results.results.PH_opt(:,chan));
            y_bl = ones(size(orig_results.results.PH_opt(:,chan))).* sin(orig_results.results.PH_opt(:,chan));
            for ct = 1:size(orig_results.results.PH_opt(:,chan),1)
                hold on; plot([0 x_bl(ct)],[0 y_bl(ct)],'ro-');
            end
            
            
            
            figPoint_ve_opt = figure;
            stem(orig_results.results.corr_mat,'filled','r.','LineWidth',1,'MarkerSize',25);
            xlabel('channels');
            ylabel('VE');
            legend('phase ref Amp','Amp','Location','NorthEast');
            title('VE opt - Amp vs ph ref Amp (UNSORTED)');
            %saveas(figPoint_ve_opt,'VE opt - Amp vs ph ref Amp','jpg');
            
            % If we want to compare data based and the model based fits
            % Load the second set of results
            if compare_phfits == 1
                [orig_fname_2, orig_fpath_2] = uigetfile('*.mat','Select original results to compare');
                orig_results_2 = load(fullfile(orig_fpath_2, orig_fname_2));
                hold on; stem(orig_results_2.results.corr_mat,'filled','k.','LineWidth',1,'MarkerSize',25);
                legend(orig_results.model.params.phase_fit,orig_results_2.model.params.phase_fit,'Location','NorthEast');
                hold off;
                
                fh_sl_diff_map=figure;
                megPlotMap(orig_results.results.corr_mat - orig_results_2.results.corr_mat,[-0.5 0.5],fh_sl_map,'jet',...
                    'VE diff',[],[],'interpmethod',interpmethod);
                
            end
                        
            [val, idx] = sort(orig_results.results.corr_mat, 'descend');
            nanChan = ~isnan(val);
            allCh = idx(nanChan);
            ch = allCh(1:15);
            timeAxis = linspace(0, length(orig_results.results.orig_model.fit_data(:,1))*1.100, length(orig_results.results.orig_model.fit_data(:,1)));
            for i =1:length(ch)
                ch_cur = ch(i);
                figPoint_1{i} = figure; set(gcf, 'Color', 'w', 'Position', [102   999   680   324])
                plot(timeAxis,orig_results.results.orig_model.fit_data(:,ch_cur),'ko-','LineWidth',1);
                hold on; plot(timeAxis,orig_results.results.orig_model.preds(:,ch_cur),'r','LineWidth',2);
                grid on; box off;
                set(gca, 'XLim', [0 timeAxis(end)], 'YLim', [-50,50].*10^-15);
                set(gca, 'YTick', linspace(-40*10^-15,40.*10^-15,5))
                xlabel('Time (s)', 'FontSize', 15)
                ylabel('Magnetic flux (Tesla)', 'FontSize', 15)
                set(gca, 'Fontsize', 15, 'TickDir', 'out');
                title(strcat('CH ',':',num2str(ch_cur),{' '},'VE',':',num2str(orig_results.results.corr_mat(ch_cur))));
            end
            
        case 'scrambled'
            
            %% for scrambled
            
            % Select scrambloing results
            [scr_fname, scr_fpath] = uigetfile('*.mat','Select scrambling results');
            scr_results = load(fullfile(scr_fpath, scr_fname));
            
            
            % Original correlation
            fh_sl_map=figure;
            megPlotMap(scr_results.results.orig_corr(:,1),[0 0.6],fh_sl_map,'jet',...
                'Phase ref fit stimulus locked',[],[],'interpmethod',interpmethod);
            %saveas(fh_sl_map,strcat('sl_map','.jpg'));
            
            % Scrambled correlation
            fh_sl_scr_map=figure;
            megPlotMap(nanmedian(scr_results.results.scrambled_corr,2),[0 0.6],fh_sl_scr_map,'jet',...
                'Phase ref fit stimulus locked scrambled',[],[],'interpmethod',interpmethod);
            %saveas(fh_sl_scr_map,strcat('sl_scr_map','.jpg'));
            
            % p value computation
            
            % Averaging channels
            
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
                
                rel = find(med_corr>0.1);
                org = find(scr_results.results.orig_corr>0.1);
                org_rel_ov = find((med_corr>0.1) == (scr_results.results.orig_corr>0.1) & (med_corr>0.1)==1);
                org_rel_noov = find((med_corr>0.1) ~= (scr_results.results.orig_corr>0.1) & (med_corr>0.1)==1);
                figPoint_rlchve_map = mprfPlotHeadLayout(good_chan);
            end
            
            n_iter = size(scr_results.results.scrambled_corr,2);
            
            VE_orig = nanmean(scr_results.results.orig_corr(good_chan,1),1);
            figPoint_scr=figure;
            scr_dist = nanmean(scr_results.results.scrambled_corr(good_chan,:),1);
            [N,x] = hist(scr_dist,100);
            hist(scr_dist,100);
            tmp_num = sum(N(find(x>VE_orig)));
            tmp_denom = size(scr_results.results.scrambled_corr,2);
            p_val = tmp_num / tmp_denom;
            xlabel('Variance explained');
            ylabel('number of outcomes');
            if p_val==0
                title(strcat('Variance explained of 1000 scrambling iteration (p < 0.001)'));
            else
                title(strcat('Variance explained 1000 scrambling iteration (p= ',num2str(p_val),')'));
            end
            hold on; plot(VE_orig*ones(1,max(N)+1),0:max(N),'r','LineWidth',2);
            %saveas(figPoint_scr,strcat('sl_scr_dist','.jpg'));
            
        case 'scrambled individual channel'
            [scr_fname, scr_fpath] = uigetfile('*.mat','Select scrambling results');
            scr_results = load(fullfile(scr_fpath, scr_fname));
            
            % Individual channels
            n_chan = size(scr_results.results.orig_corr,1);
            p_val_allchan = nan(1,n_chan);
            for cur_chan = 1:n_chan
                %                if ~isnan(nansum(scr_results.results.scrambled_corr(:,cur_chan)))
                VE_orig = scr_results.results.orig_corr(cur_chan,:);
                if ~isnan(VE_orig)
                    %figPoint_scr=figure;
                    scr_dist = scr_results.results.scrambled_corr(cur_chan,:);
                    [N,x] = hist(scr_dist,100);
                    %hist(scr_dist,100);
                    tmp_num = sum(N(find(x>VE_orig)));
                    tmp_denom = size(scr_results.results.scrambled_corr,2);
                    p_val = tmp_num / tmp_denom;
                    %xlabel('Variance explained');
                    %ylabel('number of outcomes');
                    %if p_val==0
                    %    title(strcat('Variance explained of 1000 scrambling iteration (p < 0.001)'));
                    %else
                    %    title(strcat('Variance explained 1000 scrambling iteration (p= ',num2str(p_val),')'));
                    %end
                    %hold on; plot(VE_orig*ones(1,max(N)+1),0:max(N),'r','LineWidth',2);
                    if p_val==0
                        p_val = 0.0001;
                    end
                    p_val_allchan(cur_chan) = -log10(p_val);
                end
            end
            fh_sl_scr_indchan_map=figure;
            megPlotMap(p_val_allchan,[0 5],fh_sl_scr_indchan_map,'jet',...
                'Phase ref fit stimulus locked scrambled -log10(p)',[],[],'interpmethod',interpmethod);
            
            good_chan = find(p_val_allchan > 2);
            figPoint_rlchve_map = mprfPlotHeadLayout(good_chan);
            
            
            
        case 'range'
            
            %% pRF range
            
            if strcmpi(chan_sel,'ind_scr')
                % Select scrambloing results
                [scr_fname, scr_fpath] = uigetfile('*.mat','Select scrambling results');
                scr_results = load(fullfile(scr_fpath, scr_fname));
                
                % Individual channels
                n_chan = size(scr_results.results.orig_corr,1);
                p_val_allchan = nan(1,n_chan);
                for cur_chan = 1:n_chan
                    %                if ~isnan(nansum(scr_results.results.scrambled_corr(:,cur_chan)))
                    VE_orig = scr_results.results.orig_corr(cur_chan,:);
                    if ~isnan(VE_orig)
                        %figPoint_scr=figure;
                        scr_dist = scr_results.results.scrambled_corr(cur_chan,:);
                        [N,x] = hist(scr_dist,100);
                        %hist(scr_dist,100);
                        tmp_num = sum(N(find(x>VE_orig)));
                        tmp_denom = size(scr_results.results.scrambled_corr,2);
                        p_val = tmp_num / tmp_denom;
                        %xlabel('Variance explained');
                        %ylabel('number of outcomes');
                        %if p_val==0
                        %    title(strcat('Variance explained of 1000 scrambling iteration (p < 0.001)'));
                        %else
                        %    title(strcat('Variance explained 1000 scrambling iteration (p= ',num2str(p_val),')'));
                        %end
                        %hold on; plot(VE_orig*ones(1,max(N)+1),0:max(N),'r','LineWidth',2);
                        if p_val==0
                            p_val = 0.0001;
                        end
                        p_val_allchan(cur_chan) = -log10(p_val);
                    end
                end
                
                good_chan = find(p_val_allchan > 2);
                
                
                % Select pRF range result (Size or position)
                [range_fname, range_fpath] = uigetfile('*.mat','Select pRF range results');
                range_results = load(fullfile(range_fpath, range_fname));
                par_it = range_results.model.params.sigma_range;
                
                corr_tmp=squeeze(range_results.results.corr_mat);
                n_par_it = size(corr_tmp,1);
                
                fh_sl_scr_indchan_map=figure;
                megPlotMap(p_val_allchan,[0 5],fh_sl_scr_indchan_map,'jet',...
                    'Phase ref fit stimulus locked scrambled -log10(p)',[],[],'interpmethod',interpmethod);
                %saveas(fh_sl_scr_indchan_map,strcat('sl_scr_map','.jpg'));
                
                % channles with p<0.01
                figPoint_gdch_map = mprfPlotHeadLayout(good_chan);
                %saveas(figPoint_gdch_map,strcat('sl_gdch_map','.jpg'));
                
                % graph for variance explained vs pRF size
                tmp_corr_sl = nan(n_par_it,length(good_chan));
                all_corr = range_results.results.corr_mat;
                for i=1:n_par_it
                    tmp_corr_sl(i,:) = squeeze(all_corr(:,i,good_chan,:,1));
                end
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
                
                %                 errorbar(par_it,corr_avg_sl,corr_CI_sl);
                set(gca,'TickDir', 'out');
                if strcmpi(range_results.model.type,'prf size range')
                    xlabel('Size scaling factor (ratio)');
                    set(gca,'XTick', [0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10],'XTickLabel',{'0.2', '0.5', '1.0', '1.5', '3.0', '5.0', '10'});
                    set(gca,'XScale', 'log');
                elseif strcmpi(range_results.model.type,'position (x,y) range')
                    xlabel('position (deg)');
                    set(gca,'XTick', par_it,'XTickLabel',rad2deg(range_results.model.params.x0_range));
                end
                xlim([par_it(1),par_it(end)]);
                title(sprintf('Variance explained %s sl locked', range_results.model.type));
                set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
                ylabel('Variance explained (%)');
                F = getframe(fh_sl_VEparams);
                
            else
                % Select pRF range result (Size or position)
                [range_fname, range_fpath] = uigetfile('*.mat','Select pRF range results');
                range_results = load(fullfile(range_fpath, range_fname));
                
                if strcmpi(range_results.model.type,'pRF size range')
                    par_it = range_results.model.params.sigma_range;
                elseif strcmpi(range_results.model.type,'position (x,y) range')
                    par_it = range_results.model.params.x0_range;
                end
                
                corr_tmp=squeeze(range_results.results.corr_mat);
                n_par_it = size(corr_tmp,1);
                % for tmp_cnt=1:n_par_it
                %     fh_sl_map=figure;
                %     megPlotMap(corr_tmp(tmp_cnt,:),[0 0.6],fh_sl_map,'jet',...
                %         'Phase ref fit stimulus locked',[],[],'interpmethod',interpmethod);
                %
                %     saveas(fh_sl_map,strcat('sl_map_',num2str(tmp_cnt),'.jpg'));
                % end
                
                % figPoint = openfig('Stimulus_locked_VEvsPos_occ');
                % set(gca,'XtickLabel',rad2deg(0:pi/4:2*pi))
                % saveas(figPoint,strcat('VE_Pos','.jpg'));
                
                % figPoint = openfig('Stimulus_locked_VEvsSig_occ');
                % ylim([0 0.12]);
                % saveas(figPoint,strcat('VE_Sig','.jpg'));
                
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
                    
                end
                
                % graph for variance explained vs pRF size
                tmp_corr_sl = nan(n_par_it,length(good_chan));
                all_corr = range_results.results.corr_mat;
                for i=1:n_par_it
                    tmp_corr_sl(i,:) = squeeze(all_corr(:,i,good_chan,:,1));
                end
                corr_avg_sl = nanmean(tmp_corr_sl,2);
                corr_std_sl = nanstd(tmp_corr_sl,0,2);
                corr_ste_sl = corr_std_sl ./ sqrt(size(tmp_corr_sl,2));
                corr_CI_sl = 1.96 .* corr_ste_sl;
                
                fh_sl_VEparams = figure; set(gcf, 'Position', [1000, 592, 838, 746]);
                
                % Plot mean with shaded error bar using 'patch' function
                lo = corr_avg_sl - corr_CI_sl;
                hi = corr_avg_sl + corr_CI_sl;
                color = [0.5 0.5 0.5];
                err= patch([par_it, fliplr(par_it)], [lo', fliplr(hi')], color, 'FaceAlpha', 0.5, 'LineStyle',':');               
                hold on;
                plot(par_it,corr_avg_sl,'r','Linewidth',3);
                
                % Add labels and make pretty
                set(gca,'TickDir', 'out');
                if strcmpi(range_results.model.type,'prf size range')
                    xlabel('Size scaling factor (ratio)');
                    set(gca,'XTick', [0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10],'XTickLabel',{'0.2', '0.5', '1.0', '1.5', '3.0', '5.0', '10'});
                    set(gca, 'XScale', 'log');
                elseif strcmpi(range_results.model.type,'position (x,y) range')
                    xlabel('Position (deg)');
                    set(gca,'XTick', par_it,'XTickLabel',rad2deg(range_results.model.params.x0_range));
                end
                xlim([par_it(1),par_it(end)]);
                title(sprintf('Variance explained %s sl locked', range_results.model.type));
                set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
                ylabel('Variance explained (%)');
                F = getframe(fh_sl_VEparams);

                % Plot individual sensors and their topography per
                % iteration of parameter change
                if strcmpi(range_results.model.type,'prf size range')
                    cols = 4; else, cols =2; end
                rows = round(n_par_it/cols);
                fhMeshPerParamIter = figure; set(gcf, 'Position', [326,584,1234,754], 'Color', 'w'); hold all; 
                for ii=1:n_par_it
                    mesh_data_to_plot = squeeze(all_corr(1,ii,:)); 
                    subplot(rows, cols, ii);
                    megPlotMap(mesh_data_to_plot,[0 0.6],fhMeshPerParamIter,'parula',...
                    par_it(ii),[],[],'interpmethod',interpmethod);                    
                end
                
                % Sort peaks of individual sensors
                [val_peaks, idx_peaks]  = max(tmp_corr_sl);
                [~,sensorsVEabove50prctle] = find(val_peaks>prctile(val_peaks,50));
                sensorsVEbelow50prctle = setdiff(1:length(idx_peaks), sensorsVEabove50prctle);
                
                % Plot below and above 15% variance explained in red and
                % green colors 
                fh_Ind_SL_VEparams = figure; set(gcf, 'Position', [326,584,1234,754], 'Color', 'w'); hold all;
                plot(par_it, tmp_corr_sl(:,sensorsVEabove50prctle), 'g');
                scatter(par_it(idx_peaks(sensorsVEabove50prctle)), val_peaks(sensorsVEabove50prctle), 50, 'g', 'filled');
                plot(par_it, tmp_corr_sl(:,sensorsVEbelow50prctle), 'r'); 
                scatter(par_it(idx_peaks(sensorsVEbelow50prctle)), val_peaks(sensorsVEbelow50prctle), 50, 'r', 'filled');

                set(gca,'TickDir', 'out');
                if strcmpi(range_results.model.type,'prf size range')
                    xlabel('Size scaling factor (ratio)');
                    set(gca,'XTick', [0.2, 0.5, 1.0, 1.5, 3.0, 5.0, 10],'XTickLabel',{'0.2', '0.5', '1.0', '1.5', '3.0', '5.0', '10'});
                    set(gca, 'XScale', 'log');
                elseif strcmpi(range_results.model.type,'position (x,y) range')
                    xlabel('position (deg)');
                    set(gca,'XTick', par_it,'XTickLabel',rad2deg(range_results.model.params.x0_range));
                end
                xlim([par_it(1),par_it(end)]);
                title(sprintf('Variance explained %s sl locked for individual sensors', range_results.model.type));
                set(gca, 'XGrid', 'on', 'YGrid', 'on', 'FontSize', 20); axis square;
                ylabel('Variance explained (%)');
                           
                fh_above50 = mprfPlotHeadLayout(good_chan(sensorsVEabove50prctle)', true, [], false); title('Sensors with VE above 50th percentile');
                fh_below50 = mprfPlotHeadLayout(good_chan(sensorsVEbelow50prctle)', true, [], false); title('Sensors with VE below 50th percentile');
                
            end
            
        case 'reliability'
            
            %% Reliability plot
            [rel_fname, rel_fpath] = uigetfile('*.mat','Select reliability run results');
            rel_results = load(fullfile(rel_fpath, rel_fname));
            cd(rel_fpath);
            med_corr = (nanmedian(rel_results.results.split_half.corr_mat,2));
            good_chan = find(med_corr>0.2);
            figPoint_rlve_map = mprfPlotHeadLayout(good_chan);
            %saveas(figPoint_rlve_map,strcat('rl_phrefAmp_rel_chan','.jpg'));
            
            fh_sl_rl_map=figure;
            megPlotMap(med_corr,[0 0.6],fh_sl_rl_map,'jet',...
                'reliability ph ref amplitude',[],[],'interpmethod',interpmethod);
            %saveas(fh_sl_rl_map,strcat('rl_phrefAmp_rel_map','.jpg'));
            
        case 'phase_diff'
            % Select Original results with data_fit
            [datafit_fname, datafit_fpath] = uigetfile('*.mat','Select data fit results');
            datafit_results = load(fullfile(datafit_fpath, datafit_fname));
            
            % Select Original results with model_fit
            [modelfit_fname, modelfit_fpath] = uigetfile('*.mat','Select model fit results');
            modelfit_results = load(fullfile(modelfit_fpath, modelfit_fname));
            
            datafit_PH_opt = datafit_results.results.PH_opt;
            modelfit_PH_opt = modelfit_results.results.PH_opt;
            
            diff_PH_opt = abs(wrapToPi(datafit_PH_opt - modelfit_PH_opt));
            fh_phd_map=figure;
            megPlotMap(diff_PH_opt,[0 pi],fh_phd_map,'jet',...
                'Phase diff: model and data fit',[],[],'interpmethod',interpmethod);
            
    end
    
    %% Save the figures
    
    if strcmpi(model_type,'original model')
        
        rel_dir = 'original_model';
        phase_fit = orig_results.model.params.phase_fit;
        
    elseif strcmpi(model_type,'fix prf size')
        
        rel_dir = ['fixed_prf_size_' num2str(prfsz_results.model.params.sigma_fix)];
        phase_fit = fixed_results.model.params.phase_fit;
        
    elseif strcmpi(model_type,'scrambled')
        
        rel_dir = 'scrambled';
        phase_fit = scr_results.model.params.phase_fit;
        
    elseif strcmpi(model_type,'scrambled individual channel')
        
        rel_dir = 'scrambled';
        phase_fit = scr_results.model.params.phase_fit;
        
    elseif strcmpi(model_type,'range')
        
        if strcmpi(range_results.model.type,'prf size range')
            rel_dir = 'prf_size_range';
            phase_fit = range_results.model.params.phase_fit;
            
        elseif strcmpi(range_results.model.type,'position (x,y) range')
            rel_dir = 'prf_position_range';
            phase_fit = range_results.model.params.phase_fit;
            
        end
        
    elseif strcmpi(model_type,'reliability')
        rel_dir = 'reliability';
        phase_fit = rel_results.model.params.phase_fit;
        
    elseif strcmpi(model_type,'phase_diff')
        rel_dir = 'phase';
        phase_fit = '';
        
    end
    
    type='sl';
    
    save_dir = fullfile(main_dir, res_dir, rel_dir, ['Run_' type '_' chan_sel '_' phase_fit '_' cur_time]);
    mkdir(save_dir);
    
    
    cd(save_dir);
    
    if strcmpi(model_type,'original model')
        
        saveas(fh_sl_map,strcat('sl_map','.jpg')); % original map
        
        saveas(fh_sl_ph_map,strcat('ph_map','.jpg')); % Optimum phase map
        
        saveas(figPoint_ve_opt,'VE opt - Amp vs ph ref Amp','jpg'); % VE plot
        
        for i=1:length(ch)
            ch_cur = ch(i);
            saveas(figPoint_1{i},strcat(num2str(ch_cur),'phrefAmp'),'jpg'); % example time series
        end
        
        
    elseif strcmpi(model_type,'scrambled')
        
        saveas(fh_sl_map,strcat('sl_map','.jpg')); % Original map
        
        saveas(fh_sl_scr_map,strcat('sl_scr_map','.jpg')); % Scrambled correlation map
        
        saveas(figPoint_scr,strcat('sl_scr_dist','.jpg')); % Scrambled distribution
        
        if strcmpi(chan_sel,'rel')
            saveas(figPoint_rlchve_map,strcat('sl_rel_chan','.jpg'));
        end
        
    elseif strcmpi(model_type,'scrambled individual channel')
        saveas(fh_sl_scr_indchan_map,strcat('sl_scr_indchan_p_map','.jpg')); % p value for individual channel
        
        saveas(figPoint_rlchve_map,strcat('sl_scr_indchan_goodch_map','.jpg')); % channels with p < 0.05
        
        
    elseif strcmpi(model_type,'range')
        if strcmpi(chan_sel,'ind_scr')
            
            saveas(figPoint_gdch_map,strcat('sl_gdch_map','.jpg')); % Good channels, channels with p<0.01
            
            saveas(fh_sl_scr_indchan_map,strcat('sl_scr_indchan','.jpg')); % p value
            
            imwrite(F.cdata, 'sl_VEparams_map.jpg' , 'jpg')
            
        else
            
            imwrite(F.cdata, 'sl_VEparams_map.jpg' , 'jpg');% VE vs range param (either size or position)
            saveas(fhMeshPerParamIter,strcat('sl_VEparams_individualMeshes','.jpg'));
            saveas(fh_Ind_SL_VEparams,strcat('sl_VEparams_individualSensors','.jpg'));
            saveas(fh_above50, strcat('sl_VEparams_individualSensors_above15prct','.jpg'));
            saveas(fh_below50, strcat('sl_VEparams_individualSensors_below15prct','.jpg'));

        end
        
        
        
    elseif strcmpi(model_type,'reliability')
        
        saveas(figPoint_rlve_map,strcat('rl_phrefAmp_rel_chan','.jpg'));
        
        saveas(fh_sl_rl_map,strcat('rl_phrefAmp_rel_map','.jpg'));
    end
    
    
    cd(main_dir);
    
catch ME
    if ~exist('plot_type','var') || isempty(plot_type)
        fprintf('Select one of the following options Eg: mprf_plot_figures(1) for Original model: \n')
        fprintf('1 -  Original model \n')
        fprintf('2 -  scrambled; channels for averaging  - all the channels on the back of head \n')
        fprintf('3 -  scrambled; channels for averaging based on the reliability metric \n')
        fprintf('4 -  scrambled individual channels - p value for individual channels, to select the channels that are sensitive to scrambling \n')
        fprintf('5 -  range; channels for averaging  - all the channels on the back of head \n')
        fprintf('6 -  range; channels for averaging based on the reliability metric \n')
        fprintf('7 -  range; channels for averaging are selected based on whether they are senstive to scrambling of pRF parameters \n')
        fprintf('8 -  reliability; - across repeat reliability metric \n')
        fprintf('9 -  phase_diff; - Reliable phase value for the phase referenced amplitude fits \n')
    else
        rethrow(ME)
    end
end

end