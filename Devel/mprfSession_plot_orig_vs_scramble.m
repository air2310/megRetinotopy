%mprfSession_plot_orig_vs_scramble

% Results reliability check
[rel_fname, rel_fpath] = uigetfile('*.mat','Select reliability run results');
rel_results = load(fullfile(rel_fpath, rel_fname));

rel_mat = rel_results.results.split_half.corr_mat;
med_rel = nanmedian(rel_mat,2);
good_channels = med_rel > 0.1;

% Original, bootstrapped
[orig_btstrp_fname, orig_btstrp_fpath] = uigetfile('.mat','Please select bootstrapped original model');
orig_bt_results = load(fullfile(orig_btstrp_fpath, orig_btstrp_fname));


% Original model and  Scrambled model:
[scr_fname, scr_fpath] = uigetfile('.mat','Please select scrambled run results');
scr_results = load(fullfile(scr_fpath, scr_fname));

orig_mat = scr_results.results.orig_corr;
scr_mat = scr_results.results.scrambled_corr;

gchn_orig = orig_mat(good_channels);
gchn_scr = scr_mat(good_channels,:);
gchn_scr_med = nanmedian(gchn_scr,2);
gchn_scr_low = prctile(gchn_scr,2.5,2);
gchn_scr_high = prctile(gchn_scr,97.5,2);
gchn_orig_bt = orig_bt_results.results.corr_ci_sl(:,good_channels)';

[~, sort_idx] = sort(gchn_orig);

x_axes_labels = find(good_channels);
x_axes = 1:sum(good_channels);

figure; hold on;
p1 = plot(x_axes, gchn_orig(sort_idx),'k-','LineWidth',2);
p2 = plot(x_axes, gchn_scr_med(sort_idx),'r-');
p3 = plot(x_axes, [gchn_scr_low(sort_idx) gchn_scr_high(sort_idx)],'r--');
set(gca,'XTick',x_axes,'XTickLabel',x_axes_labels(sort_idx))
grid on
legend([p1 p2],'Original model','Scrambled model','Location','NW');
xlabel('Channel number')
ylabel('Variance explained')

[~, sort_idx] = sort(gchn_orig_bt(:,2));

figure; hold on; 
p1 = plot(x_axes, gchn_orig_bt(sort_idx,2),'k-','LineWidth',2);
p2 = plot(x_axes, gchn_orig_bt(sort_idx,[1 3]),'k--','LineWidth',1);
p3 = plot(x_axes, gchn_scr_med(sort_idx),'r-');
set(gca,'XTick',x_axes,'XTickLabel',x_axes_labels(sort_idx))
grid on
legend([p1 p3],'Original model','Scrambled model','Location','NW');
xlabel('Channel number')
ylabel('Variance explained')











figure; hold on;
plot(median(rel_mat,2))
plot(orig_mat,'r--')




return
scan_mat = rel_results.results.scans.corr_mat;

scan_mat2 = nan(157,18);

for n = 2:19
    for nn = 1:157
   scan_mat2(nn,n-1) = median(scan_mat{n}(nn,:)); 
    
    end
    
end

scan_mat3 = squeeze(nanmedian(scan_mat2,2));
figure; plot(scan_mat2(103,:))












% Position of the sensors
mprfPlotHeadLayout;