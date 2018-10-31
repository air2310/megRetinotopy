function [] = loadAndVisualizeStim(stimPths, stimType)


switch stimType
    
    case 'MEG'
        
        % Load stimuli and grid
        stim = load(stimPths.MEGStim.pth);
        stim = stim.stimulus;
        
        grid = load(stimPths.MEGStimGrid.pth);
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
        
        % [EK]: Why is [101 by 101] the image output size?
        imOutSize = [101 101];
        megStim.resizedIm = zeros([imOutSize, length(stimTriggers)]);
        
        % [EK]: Why re-sample grid with rate = .2105???
        newGridSampleRate = max(grid.Xd(:)) ./ floor(max(imOutSize) / 2);
        
        % Loop over triggers,
        for n = 1:length(stimTriggers)
            % Get stimulus images
            currIm = double(stim.images(:,:,stimTriggers(n)));
            
            % Resize stimulus to new size
            currIm = imresize(currIm,imOutSize,'nearest');
            
            % Convert to resized binary mask
            megStim.resizedIm(:,:,n) = double(ceil(abs(currIm - cGray) ./ max(abs(cRange - cGray)))).*(newGridSampleRate.^2);
            
        end
        
        % Get new resized X and Y dimensions
        [megStim.resizedX, megStim.resizedY] = meshgrid(linspace(xRange(1), xRange(2), imOutSize(1)),...
            linspace(yRange(1), yRange(2), imOutSize(2)));
        
        % Get full field binary mask, given all stimulus positions
        fullFieldBinaryMask = sum(megStim.resizedIm,3) > 0;
        
        % Concatenate first two dimensions of image (101 x 101 --> 10201),
        % by number of stimulus frames (211)
        megStim.im = reshape(megStim.resizedIm,[],size(megStim.resizedIm,3));
        
        % Mask images with binary mask --> [EK]: why? Wasn't it already
        % binary, and aren't they the same size?
        megStim.im = megStim.im(fullFieldBinaryMask(:),:);
        
        % Add concatenated image x and y size to megStim struct
        megStim.X  = megStim.resizedX(:);
        megStim.Y  = megStim.resizedY(:);
        megStim.X  = megStim.X(fullFieldBinaryMask);
        megStim.Y  = megStim.Y(fullFieldBinaryMask);
        
        
%             sz1 = numel(X);
%             sz2 = 
%             X   = repmat(X(:),1,sz2);
%             Y   = repmat(Y(:),1,sz2);
% 
%             define centered around (0,0) coordinate space
%             x0 = repmat(x0(:)',sz1,1);
%             y0 = repmat(y0(:)',sz1,1);
% 
%             Translate grid so that center is at RF center
%             X = X - x0;   % positive x0 moves center right
%             Y = Y - y0;   % positive y0 moves center up
        
    case 'MRI'
        
        prfParams = load(stimPths.PRFParams.pth);
        
        rm_stim.im = prfParams.stim.images_unconvolved;
        rm_stim.im_conv = prfParams.analysis.allstimimages';
        
        rm_stim.window = prfParams.stim.stimwindow;
        rm_stim.X = prfParams.analysis.X;
        rm_stim.Y = prfParams.analysis.Y;
        
        [rm_stim.full_x, rm_stim.full_y] = meshgrid(unique(rm_stim.X), unique(rm_stim.Y));
        [~, rm_stim.full_im] = rmStimulusMatrix(rmparams,[],[],false,false);
        
        
end