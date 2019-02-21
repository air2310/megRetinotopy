% s_visualizeROIvsVEmask

subjectSessionDir = '/Volumes/server/Projects/MEG/Retinotopy/Subject_sessions';
bsDir = '/Volumes/server/Projects/MEG/brainstorm_db';
subject = 'wlsubj068';
anatDir = fullfile(bsDir, 'MEG_Retinopy', 'anat', subject);

cd(fullfile(subjectSessionDir, subject, 'mask'))


tmp = load('roimask.mat');
m.roi = double(tmp.cur_in);
tmp = load('vemask_withroimask.mat');
m.both = double(tmp.cur_in);
tmp = load('vemask_noroimask.mat');
m.ve = double(tmp.cur_in);

clear tmp;



visualizeBrainstormMesh(anatDir, m.both, [],[],[], 'ROI+VE mask');
print(gcf,'-dpng', [] ,fullfile(subjectSessionDir, subject, 'mask', 'bs_mesh_roive'));
savefig(gcf,fullfile(subjectSessionDir, subject, 'mask', 'bs_mesh_roive.fig'));

visualizeBrainstormMesh(anatDir, m.ve, [],[],[], 'VE mask');
print(gcf,'-dpng', [] ,fullfile(subjectSessionDir, subject, 'mask', 'bs_mesh_ve'));
savefig(gcf,fullfile(subjectSessionDir, subject, 'mask', 'bs_mesh_ve.fig'));


fprintf('Number of total vertices ROI+VE mask: %d\n', sum(m.both));
fprintf('Number of total vertices VE mask: %d\n', sum(m.ve));
fprintf('Number of total vertices ROI mask: %d\n', sum(m.roi));
fprintf('Number of vertices overlapping ROI+VE vs VE only mask: %d\n', length(intersect(find(m.both),find(m.ve))));
fprintf('Number of vertices overlapping ROI+VE vs ROI only mask: %d\n', length(intersect(find(m.both),find(m.roi))));
fprintf('Number of vertices overlapping ROI only vs VE only mask: %d\n', length(intersect(find(m.roi),find(m.ve))));