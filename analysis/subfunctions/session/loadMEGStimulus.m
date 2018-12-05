function meg_stim = loadMEGStimulus(s)

% subfunction of loadStim, to specifically prepare the MEG stimulus

% Load stimuli and grid
stim = load(s.MEGStim.pth);
stim = stim.stimulus;

grid = load(s.MEGStimGrid.pth);
grid = grid.grid;

% Get x and y dim range from grid
xRange = [min(grid.Xd(:)) max(grid.Xd(:))];
yRange = [min(grid.Yd(:)) max(grid.Yd(:))];

% Get min and max image values
cRange = double([min(stim.images(:)) max(stim.images(:))]);

% Get mid value of RGB for blank (128)
im1 = stim.images(:,:,1);
cGray = double(mode(im1(:)));

% Get triggers for stim and corresponding condition number
stimTriggers = stim.seq(stim.trigSeq > 0);

% We use [101 by 101] as the image output size, because the
% prf parameters has a number of stimulus grid points of 50,
% and length(-50:50)= 101. So the X and Y grid used when solving
% the retinotopy model with MRI data was 101x101. See vistasoft
% function: makeStimFromScan(params,1);
prfParams      = load(s.PRFParams.pth);
szDownSampled  = sqrt(size(prfParams.params.stim.stimwindow,1));
imOutSize      = [szDownSampled szDownSampled]; % imOutSize = [101 101];

meg_stim.resizedIm = zeros([imOutSize, length(stimTriggers)]);

% Redefine the grid with sample rate with the used visual angle
newGridSampleRate = max(grid.Xd(:)) ./ floor(max(imOutSize) / 2);

% Loop over triggers,
for n = 1:length(stimTriggers)
    % Get stimulus images
    currIm = double(stim.images(:,:,stimTriggers(n)));
    
    % Resize stimulus to new size
    currIm = imresize(currIm,imOutSize,'nearest');
    
    % Convert to resized binary mask
    meg_stim.resizedIm(:,:,n) = double(ceil(abs(currIm - cGray) ./ max(abs(cRange - cGray)))).*(newGridSampleRate.^2);
    
end

% Get new resized X and Y dimensions
[meg_stim.resizedX, meg_stim.resizedY] = meshgrid(linspace(xRange(1), xRange(2), imOutSize(1)),...
    linspace(yRange(1), yRange(2), imOutSize(2)));

% Get full field binary mask, given all stimulus positions
meg_stim.window = sum(meg_stim.resizedIm,3) > 0;

% Concatenate first two dimensions of image (101 x 101 --> 10201),
% by number of stimulus frames (211)
meg_stim.im = reshape(meg_stim.resizedIm,[],size(meg_stim.resizedIm,3));

% Mask images with binary mask --> [EK]: why? Wasn't it already
% binary, and aren't they the same size? (Note: previously called full_im)
meg_stim.im = meg_stim.im(meg_stim.window(:),:);

% Add concatenated image x and y size to megStim struct
meg_stim.X  = meg_stim.resizedX(:);
meg_stim.Y  = meg_stim.resizedY(:);
meg_stim.X  = meg_stim.X(meg_stim.window);
meg_stim.Y  = meg_stim.Y(meg_stim.window);

saveDir = fullfile(s.outPut.pth, 'stim');
if ~exist('saveDir', 'dir'); mkdir(saveDir); end;
save(fullfile(saveDir, 'meg_stimulus'),'meg_stim')

end