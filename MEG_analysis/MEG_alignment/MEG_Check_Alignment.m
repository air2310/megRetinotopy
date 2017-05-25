marker_file = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wl_subj_004/Raw/Export_02/R0774_Marker1_5.12.17.sqd';
HS_file = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wl_subj_004/Raw/Export_02/R0774_5.12.17_HS.txt';
Points_file = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wl_subj_004/Raw/Export_02/R0774_5.12.17_Points.txt';

labels = {'LPA', 'RPA', 'CPF', 'LPF', 'RPF'};

fid = fopen(Points_file, 'r');
txtCell = textscan(fid, '%f%f%f', 'Delimiter', '\n', 'CommentStyle', '%');
fclose(fid);
points_xyz = [txtCell{1}, txtCell{2}, txtCell{3}]' ./ 1000;
points_xyz = points_xyz(:,4:end);

fid = fopen(HS_file, 'r');
txtCell = textscan(fid, '%f%f%f', 'Delimiter', '\n', 'CommentStyle', '%');
fclose(fid);
hs_xyz = [txtCell{1}, txtCell{2}, txtCell{3}]' ./ 1000;

step_size = round(size(hs_xyz,2) / 1000);
hs_xyz = hs_xyz(:,1:step_size:end);

header.coreg = getYkgwHdrCoregist(marker_file);
xyzMeg = cat(1, header.coreg.hpi.meg_pos)';



[origR,origT] = rot3dfit(points_xyz', xyzMeg');
origR = origR';
origT = origT';

orig_Tmat = [origR, origT; 0 0 0 1];




%%% POSSIBLE FIX:
scale_n  = 10^-2;
% Use affine fit to find the normal of the plane through the head points:
[n_orig, v, p] = affine_fit(points_xyz');

% Normal is normalized, so rescale:
n = n_orig * scale_n;

% Plot the head points:
figure; 
plot3(points_xyz(1,:) ,points_xyz(2,:), points_xyz(3,:),...
    'R.','MarkerSize',24)
hold on

% Plot the head shape:
plot3(hs_xyz(1,:),hs_xyz(2,:), hs_xyz(3,:), 'k.');

% Compute the mean position of the head shape:
m_hs = mean(hs_xyz,2);


% Determine on which side of the plane the head shape is by finding the
% normal that closest the the head shape's mean:
dist_01 = sum(((p + -n') - m_hs').^2);
dist_02 = sum(((p + n') - m_hs').^2);

% Pick the best:
if dist_01 < dist_02
    n = -n;
end

hdl_point = p + n';
quiver3(p(1),p(2),p(3),n(1),n(2),n(3),'b')
plot3(hdl_point(1),hdl_point(2),hdl_point(3),'b^','MarkerFaceColor','b','MarkerSize',12)
points_xyz_hdl = [points_xyz hdl_point'];


scale_n2  = 10^-2;

% Do the same thing for the MEG marker points:
[n2_orig, v2, p2] = affine_fit(xyzMeg');
n2 = n2_orig.*scale_n2;

% Add points along the normal orthogonal to the plane in both directions:
hndl_01 = p2 + -n2';
hndl_02 = p2 + n2';



% Plot the head points:
figure; 
plot3(xyzMeg(1,:) ,xyzMeg(2,:), xyzMeg(3,:),...
    'R.','MarkerSize',24)
hold on
quiver3(p2(1),p2(2),p2(3),n2(1),n2(2),n2(3),'r')
quiver3(p2(1),p2(2),p2(3),-n2(1),-n2(2),-n2(3),'b')

xyzMEG_hdl_01 = [xyzMeg hndl_01'];
xyzMEG_hdl_02 = [xyzMeg hndl_02'];


% Recompute the transformation:
% This one is almost the same as the initial transformation:
[R_01, T_01] = rot3dfit(points_xyz_hdl', xyzMEG_hdl_01');
R_01 = R_01';
T_01 = T_01';

Tmat_01 = [R_01, T_01; 0 0 0 1];

% This one is not, and when replace with the transformation matrix
% Brainstorm computes, has a much better result:
[R_02, T_02] = rot3dfit(points_xyz_hdl', xyzMEG_hdl_02');
R_02 = R_02';
T_02 = T_02';

Tmat_02 = [R_02, T_02; 0 0 0 1];


t_points_01 = Tmat_01(1:3,1:3) * points_xyz;
t_points_01 = bsxfun(@plus, t_points_01, Tmat_01(1:3,4));

t_points_02 = Tmat_02(1:3,1:3) * points_xyz;
t_points_02 = bsxfun(@plus, t_points_02, Tmat_02(1:3,4));



figure; hold on;
plot3(t_points_01(1,:), t_points_01(2,:), t_points_01(3,:),'r.','MarkerSize',24)
t_01_ax = axis;

plot3(t_points_02(1,:), t_points_02(2,:), t_points_02(3,:),'g.','MarkerSize',24)
t_02_ax = axis;

plot3(xyzMeg(1,:), xyzMeg(2,:), xyzMeg(3,:),'k.','MarkerSize',24)
meg_ax = axis;




figure;
text(t_points_01(1,:), t_points_01(2,:), t_points_01(3,:),labels')
axis(t_01_ax)

figure;
text(t_points_02(1,:), t_points_02(2,:), t_points_02(3,:),labels')
axis(t_02_ax)

figure;
text(xyzMeg(1,:), xyzMeg(2,:), xyzMeg(3,:),labels')
axis(meg_ax)












