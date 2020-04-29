% Load all the coordinates, as imported by brainstorm:
MEG_alignment_issues_02;

% Adjust the scale of the head shape:
head_shape = head_shape ./ 1000;

% Compute the transformation between head points and MEG markers:
[R_wl,T_wl] = rot3dfit(wl_subj_40_points(4:end,:), wl_subj_40_markers);
R_wl = R_wl';
T_wl = T_wl';

% Report in 4x4 format:
wl_tmat = [R_wl, T_wl; 0 0 0 1];

% Do the same for Ernies data:
[R_ernie,T_ernie] = rot3dfit(ernie_points(4:end,:), ernie_markers);
R_ernie = R_ernie';
T_ernie = T_ernie';

% Report in 4x4 format
ernie_tmat = [R_ernie, T_ernie; 0 0 0 1];

% Use affine fit to find the normal of the plane through the head points:
[n, v, p] = affine_fit(wl_subj_40_points(4:end,:));

% Normal is normalized, so rescale:
n = n./100;

% Plot the head points:
figure; 
plot3(wl_subj_40_points(4:end,1) ,wl_subj_40_points(4:end,2), wl_subj_40_points(4:end,3),...
    'R.','MarkerSize',24)
hold on

% Plot the normal of the surface through the plane, in both directions:
quiver3(p(1),p(2),p(3),n(1),n(2),n(3),'b')
quiver3(p(1),p(2),p(3),-n(1),-n(2),-n(3),'r')

% Plot the head shape:
plot3(head_shape(:,1),head_shape(:,2), head_shape(:,3), 'k.');

% Compute the mean position of the head shape:
m_hs = mean(head_shape);

% Plot mean position of the head shape:
plot3(m_hs(1), m_hs(2), m_hs(3), 'bd','MarkerFaceColor','b','MarkerSize',18)

% Determine on which side of the plane the head shape is by finding the
% normal that closest the the head shape's mean:
dist_01 = sum(((p + -n') - m_hs).^2);
dist_02 = sum(((p + n') - m_hs).^2);

% Pick the best:
if dist_01 < dist_02
    n = -n;
end

% Add this point to the head points:
hndl = p + n';
wl_subj_40_points_hndl = [wl_subj_40_points(4:end,:);  hndl];

% Do the same thing for the MEG marker points:
[n2, v2, p2] = affine_fit(wl_subj_40_markers);
n2 = n2./100;

% Add points along the normal orthogonal to the plane in both directions:
hndl_01 = p2 + -n2';
hndl_02 = p2 + n2';

wl_subj_40_markers_hndl_01 = [wl_subj_40_markers ;hndl_01];
wl_subj_40_markers_hndl_02 = [wl_subj_40_markers ;hndl_02];

% Recompute the transformation:
% This one is almost the same as the initial transformation:
[R_01, T_01] = rot3dfit(wl_subj_40_points_hndl, wl_subj_40_markers_hndl_01);

% This one is not, and when replace with the transformation matrix
% Brainstorm computes, has a much better result:
[R_02, T_02] = rot3dfit(wl_subj_40_points_hndl, wl_subj_40_markers_hndl_02);





