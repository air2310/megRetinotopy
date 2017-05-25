
% Assumes vistadisp is on you Matlab path and you're on the meg_compatible
% branch. Note to self: use git checkout master and git checkout
% meg_compatible
% 
% addpath(genpath('~/Documents/MATLAB/toolboxes/vistadisp'));
% 
% % Needs psychtoolbox 3
% tbUse('psychtoolbox-3');



% Get screen parameters usefull for knowing degrees visual angle
cal = 'meg_lcd';
d   = loadDisplayParams(cal);

d.pixelSize = min(d.dimensions ./ d.numPixels); % Horizontal and vertical dimensions differ a lot, check this with someone...

% Define the orientations here. These are the orientations for the regular
% pRF mapping stimulus. There are only four, we duplicate, flip and add
% these to the first four:
orientations = (0:45:360)./360*(2*pi); % degrees -> rad
orientations = orientations([1 6 3 8]);

% Outer radius of the stimulus window:
outerRad                = 12;

% Frame rate of the monitor/projector. Needed for the frame timing
% sequence:
frame_rate              = 60;

% Changes per second:
flicker_freq            = 20;

% Duration of one stimulus frame:
stimframe               = (1/frame_rate) * (frame_rate/flicker_freq);

% Duration of one bar step:
steps_per_sec           = 1; 

% How many steps for each bar crossing, by which I mean, how many bar
% positions are there for one crossing. The amount of steps is this number
% -1
n_steps_per_pass        = 20;   % step_nx = duration.cycle.seconds./params.tr/8;

% How many images per complete cycle of checks, as we want a flickering
% stimulus, this needs to be 2
motionSteps             = 2;    

% Width of the bar:
bar_width               = 2; 

% Cycles per bar, i.e. how many black and white checks across the bar's width
cycles_per_bar          = 5;            

minCmapVal              = 0;            % minimum color map value -> might be 1 as well
maxCmapVal              = d.maxRgbValue;    % Maximum color map value -> might be 254 as well
bk                      = 128;      % Back ground

% Amount of images to be made:
numImages               = length(orientations) * n_steps_per_pass;

% size of the stimulus matrix in pixels, and we are using a square:
m = round(angle2pix(d, outerRad*2));
n = m;

% X, Y and eccentricity/radius for stimulus window:
[x,y]=meshgrid(linspace(-outerRad,outerRad,n),linspace(outerRad,-outerRad,m));
r = sqrt (x.^2  + y.^2);

% Determine when the bar orientation needs to change in the stimulus matrix:
remake_xy    = zeros(1,numImages)-1;
remake_xy(1:length(remake_xy)/length(orientations):length(remake_xy)) = orientations;

% Duplicate x and Y:
original_x   = x;
original_y   = y;

step_x       = (2*outerRad) ./ n_steps_per_pass;
step_startx  = (n_steps_per_pass-1)./2.*-step_x - (bar_width./2);


images=zeros(m,n,numImages*motionSteps,'uint8');

for imgNum=1:numImages
    
    if remake_xy(imgNum) >=0,
        x = original_x .* cos(remake_xy(imgNum)) - original_y .* sin(remake_xy(imgNum));
        y = original_x .* sin(remake_xy(imgNum)) + original_y .* cos(remake_xy(imgNum));
        % Calculate checkerboard.
        % Wedges alternating between -1 and 1 within stimulus window.
        % The computational contortions are to avoid sign=0 for sin zero-crossings
        wedges    = sign(round((cos((x+step_startx)*cycles_per_bar*(2*pi/bar_width)))./2+.5).*2-1);
        posWedges = find(wedges== 1);
        negWedges = find(wedges==-1);
        rings     = zeros(size(wedges));
        
        checks    = zeros(size(rings,1),size(rings,2),motionSteps);
        for ii=1:motionSteps,
            tmprings1 = sign(2*round((cos(y*cycles_per_bar*(2*pi/bar_width)+(ii-1)/motionSteps*2*pi)+1)/2)-1);
            tmprings2 = sign(2*round((cos(y*cycles_per_bar*(2*pi/bar_width)-(ii-1)/motionSteps*2*pi)+1)/2)-1);
            rings(posWedges) = tmprings1(posWedges);
            rings(negWedges) = tmprings2(negWedges);
            
            checks(:,:,ii)=minCmapVal+ceil((maxCmapVal-minCmapVal) * (wedges.*rings+1)./2);
        end;
        
        % reset starting point
        loX = step_startx - step_x;
    end;
    
    
    loX   = loX + step_x;
    hiX   = loX + bar_width;
    
    % This isn't as bad as it looks
    % Can fiddle with this to clip the edges of an expanding ring - want the ring to completely
    % disappear from view before it re-appears again in the middle.
    
    % Can we do this just be removing the second | from the window
    % expression? so...
    window = ( (x>=loX & x<=hiX) & r<outerRad);
    
    % yet another loop to be able to move the checks...
    
    tmpvar = zeros(m,n);
    tmpvar(window) = 1;
    tmpvar = repmat(tmpvar,[1 1 motionSteps]);
    window = tmpvar == 1;
    img         = bk*ones(size(checks));
    img(window) = checks(window);
    images(:,:,(imgNum-1).*motionSteps+1:imgNum.*motionSteps) = uint8(img); %#ok<*BDSCA>
    
    
    fprintf('.');drawnow;
end
fprintf('Done.\n');

images_all = cat(3, images, flip(flip(images,1),2), uint8(zeros(m,n)));



total_sec = steps_per_sec * n_steps_per_pass * 8;
n_frames = total_sec * flicker_freq;
seqtiming = (0:n_frames-1) .* stimframe;
n_im_per_step =  steps_per_sec * flicker_freq /2;

n_images = size(images_all,3)-1;
aa = mod(2:n_frames+(2-1),2)+1;

n_steps = 8 * n_steps_per_pass;

bb = repmat(0:2:(2*n_steps)-1,(n_frames / (n_images/2)),1);
cc = bb(:);
dd = aa(:) + cc;




idx1 = 0;
figure; 
for n = 1:length(dd); 
    imagesc(images_all(:,:,dd(n))); 
    colormap gray; 
    drawnow; 
    pause(1/flicker_freq); 
end









