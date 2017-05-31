function msh = mprfSessionColorMesh(msh, data, cmap, drange, mask)
%% SOME COMMENT TO CHECK IT
prefs = mrmPreferences;

windowID = meshGet(msh,'windowid');


sz = size(meshGet(msh,'colors'),2);
colors = ones(4, sz) .* 127;

colors(1:3,:) = meshData2Colors(data, cmap', drange, 1);
colors = rescale2(colors, [0 1], [0 255]);

% Assign the anatomy colors (usually representing curvature) to the
% locations where there are no data values.
anatColors = meshCurvatureColors(msh);
colors(1:3,~mask) = anatColors(1:3,~mask);

% Manually apply an alpha channel the overlay colors- adjusting them
% to allow the underlying mesh colors to come through. (Ie. simulate
% transparency). This is useful when your mesh is completely painted, but
% you need some clue about the surface curvature (which is presumable what
% the original mesh colors represent). -Bob
if prefs.overlayModulationDepth > 0
    a = prefs.overlayModulationDepth; % 'alpha' level
    colors(:,mask) = ((1-a) * colors(:,mask)) + (a * double(anatColors(:,mask)));
end

colors(colors>255) = 255;
colors(colors<0) = 0;
colors = uint8(round(colors))';

% Now we compute the color for every mesh vertex based on the assigned pRF
% parameter values:
%[mrv_clrs, mrv_msh] = meshOverlay(mrv_msh,surf_pa,[],hsv(256),'auto','showdata','phasedata');

% And visualize the mesh if it is not visualized yet:


if windowID == -1
    msh = meshVisualize(msh);
end

% Color it
msh = mrmSet(msh, 'colors', colors);

end



