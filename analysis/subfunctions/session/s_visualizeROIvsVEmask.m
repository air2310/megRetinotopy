% s_visualizeROIvsVEmask

tmp = load('roimask.mat');
m.roi = tmp.roi.mask;
tmp = load('ve_mask.mat');
m.ve = tmp.cur_in;

tmp = load('roimask_noroimask.mat');
m.all = tmp.roi.mask;
tmp = load('ve_mask_noroimask.mat');
m.ve_noroi = tmp.cur_in;

clear tmp;

bsDir = '/Volumes/server/Projects/MEG/brainstorm_db';
subject = 'wlsubj040';
anatDir = fullfile(bsDir, 'MEG_Retinopy', 'anat', subject);


visualizeBrainstormMesh(anatDir, m.roi, [],[],[], 'ROI mask');
visualizeBrainstormMesh(anatDir, double(m.ve), [],[],[], 'VE mask');


fprintf('Number of total vertices ROI mask: %d\n', sum(m.roi));
fprintf('Number of total vertices VE mask: %d\n', sum(m.ve));
fprintf('Number of vertices overlapping ROI/VE mask: %d\n', length(intersect(find(m.roi),find(m.ve))));

fprintf('Number of total vertices ROI mask: %d\n', sum(m.all));
fprintf('Number of total vertices VE mask: %d\n', sum(m.ve_noroi));
fprintf('Number of vertices overlapping ROI/VE mask: %d\n', length(intersect(find(m.all),find(m.ve_noroi))));