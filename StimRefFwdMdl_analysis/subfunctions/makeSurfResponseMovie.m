function makeSurfResponseMovie(predSurfResponse, meg_stim, dirPth, opt)

% Define image size
scr_sz = get(0,'ScreenSize');
rows = ceil(scr_sz(3)/3);
cols = ceil(scr_sz(3)/3);
rgb  = 3;
nTimePoints =  size(predSurfResponse,1);

fH = figure(1); clf; set(gcf, 'Position', [1 1 rows, cols], 'Color', 'w');

% Allocate space for movie frames
% meshim = zeros(rows, cols, rgb, nTimePoints, 'uint8');

% Set colormap, colormap limits and colormap threshold
cmap     = hot(256);
clim     = [min(predSurfResponse(:)) max(predSurfResponse(:))];
thresh   = .75;
meshType = 'smooth'; % 'smooth' for inflated, 'unsmooth' for pial

fprintf('\n(%s): Creating frames.', mfilename);

for t = 1:nTimePoints
    fprintf('.')
    tmp = visualizeBrainstormMesh(dirPth.bs.anatPth, predSurfResponse(t,:), cmap, thresh, clim, meshType, num2str(t));
    meshim(:,:,:,t) = uint8(tmp.cdata);
end

% Get size of mesh frame and use it to resize stim image
sz = size(meshim,1);
stim = imresize(meg_stim.resizedIm, [sz sz]);

% Make it a rgb image, even though its grayscale
stim = cat(4, stim, stim, stim); stim = permute(stim, [1 2 4 3]);

% Set value range from 1-255
mx = max(stim(:));
stim = uint8(stim / mx * 255);

% Concatenate stimulus and mesh frames
movieim = cat(2, stim, meshim);
 
% Create movie
mov = implay(movieim, 4);

if opt.saveFig
    fname = fullfile(dirPth.model.saveFigPth, opt.subfolder, 'prfSurfResponse');
    vid = VideoWriter(fname, 'MPEG-4');
    open(vid);
    
    writeVideo(vid, movieim);
end

fprintf('\n(%s): Done!', mfilename)