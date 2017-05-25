% MEG RETINOTOPY
%   MAKE CARRIER PATTERNS
%
% GOAL: Make multiple exemplars or several spatial patterns for testing
% their effectiveness in driving BOLD responses in visual cortex:
%
% (1) Narrow-band curvy patterns (ala SOC paper)
% (2) Dartboards
% (3) Checkerboards
%
% I have build this as an on/off stimulus to be analysed with corAnal

%% Stimulus parameters:

% General stimulus:
n_stim          = 3;    % one for every pattern
im_labels       = {'patterns','dartboard','checkerboard'}; % Names of patterns for saving
n_repeats       = 6;    % cycles per stimulus
n_fix_changes   = 20;   % amount of fixation color changes
n_runs          = 5;
period          = 12; % Seconds stim on/off
flicker_freq    = 10; % Hz
save_path       = '/Volumes/server/Projects/MEG/Retinotopy/Stimuli';
TR              = 1.5; % Seconds

if mod(period,TR)
    error('Period is not a whole multiple of TR')
    
else
    period_tr       = period / TR;
    
end

% Screen options
cal             = 'CBI_NYU_projector';
d               = loadDisplayParams(cal);
d.pixelSize     = min(d.dimensions ./ d.numPixels);
res             = min(d.numPixels); % Size of stimulus in pixels. Stimulus will seen through circular window with this diameter
frame_rate      = d.frameRate;
screen_size     = 2*pix2angle(d,res/2);

% Curvy patterns options
num_examples    = 32;   % Amount of patterns generated
sparseness      = 1/10; % Sparseness of curvy patterns (bigger denominator = more sparse)
center_sf       = 3;    % cycles per degree
bandwidth       = 2;    % Bandwidth of filter in octaves
filter_size     = 31;   % size of filter in pixels
w_filter        = 'flt_dog';

% Addpaths
addpath(genpath('~/matlab/git/soccode')) % Add paths and load the bpfilter and the base images
addpath(genpath('~/matlab/git/knkutils'))    
addpath(genpath('~/matlab/git/matlabPyrTools'))

% Dartboard options
wedges_dpc      = 2;    % i.e. degrees for one black/white repeat 
rings_dpc       = 2;

n_wedges        = screen_size / wedges_dpc * 2; % Amount of wedges making up the dartboard
n_rings         = screen_size / rings_dpc * 2;  % Amount of rings

% Checkerboard options
rows_dpc        = 2;
cols_dpc        = 2;

n_rows          = screen_size / rows_dpc * 2; % Amount of checks along the board vertically
n_cols          = screen_size / cols_dpc * 2; % Amount of checks along the board horizontally

% Store general options
opts.general.frame_rate = frame_rate;
opts.general.res = res;
opts.general.flicker_freq = flicker_freq;
opts.general.screen_size = screen_size;
opts.general.TR = TR;
opts.general.period_tr = period / TR;

% Store curvy pattern options:
opts.patterns.num_examples = num_examples;
opts.patterns.sparseness = sparseness;
opts.patterns.filter_center = center_sf;
opts.patterns.filter_bw = bandwidth;
opts.patterns.filter_size = filter_size;

% Store dartboard options:
opts.dartboard.wedges_dpc = wedges_dpc;
opts.dartboard.rings_dpc = rings_dpc;
opts.dartboard.n_wedges = n_wedges;
opts.dartboard.n_rings = n_rings;

% Store checkerboard options:
opts.checkerboard.rows_dpc = rows_dpc;
opts.checkerboard.cols_dpc = cols_dpc;
opts.checkerboard.n_rows = n_rows;
opts.checkerboard.n_cols = n_cols;


%% General stuff

[X, Y] = meshgrid(-(res-1)/2:(res-1)/2);
[TH, R] = cart2pol(X,Y);
TH(TH < 0) = TH(TH < 0) + 2 * pi;

window = R < (res)/2;

for run = 1:n_runs
    %% (1) Curvy patterns
    
 
    try
        [flt_cos, flt_dog] = mprfComputeAndStoreBPFilter('mri_filter');
        
    catch
        [flt_cos, flt_dog] = mprfComputeAndStoreBPFilter([], center_sf, bandwidth, filter_size, res, screen_size, 'mri_filter',false);

    end

    im_01 = zeros(res, res, num_examples);
    
    for ii = 1:num_examples
        im_01(:,:,ii) = createPatternStimulus([res, res], 1/10, eval(w_filter)) .* window;
    
    end
    
    im_01 = im_01 + 0.5;
    im_01(im_01 < 0) = 0;
    im_01 = uint8(im_01 .* 255);
    
    
    %% (2) Dartboard
    
    im_02 = nan(res,res,2*n_repeats);
    
    for n = 1:2:n_repeats*2
        
        rph = rand .* (2*pi);
        
        
        TH(TH < 0) = TH(TH < 0) + 2 * pi;
        wedges = sign(2*round((sin(TH .* (n_wedges/2) + rph)+1)/2)-1);
        rings_01 = sign(2*round((sin( (R./((res-1)/2)).*2*pi .* (n_rings/2) +rph )+1)/2)-1);
        rings_02 = sign(2*round((sin( (R./((res-1)/2)).*2*pi .* (n_rings/2) +pi +rph )+1)/2)-1);
        
        rings = zeros(res,res);
        rings(wedges < 0) = rings_01(wedges < 0);
        rings(wedges > 0) = rings_02(wedges > 0);
        
        im_02(:,:,n) = uint8(((rings / 2) .* window + 0.5) .* 255);
        
        rings = zeros(res,res);
        rings(wedges < 0) = rings_02(wedges < 0);
        rings(wedges > 0) = rings_01(wedges > 0);
        
        im_02(:,:,n+1) = uint8(((rings / 2) .* window + 0.5).*255);
    end
    
    
    %% (3) Checkerboard
    im_03 = nan(res,res,2*n_repeats);
    
    for n = 1:2:n_repeats*2
        rph = rand .* (2*pi);
        
        rows = sign(2*round((sin((Y./((res-1)/2)).*2*pi .* (n_rows/4) + rph)+1)/2)-1);
        cols_01 = sign(2*round((sin( (X./((res-1)/2)).*2*pi .* (n_cols/4) +rph )+1)/2)-1);
        cols_02 = sign(2*round((sin( (X./((res-1)/2)).*2*pi .* (n_cols/4) +pi +rph )+1)/2)-1);
        
        cols = zeros(res,res);
        cols(rows < 0) = cols_01(rows < 0);
        cols(rows > 0) = cols_02(rows > 0);
        
        im_03(:,:,n) = uint8(((cols / 2) .* window + 0.5) .* 255);
        
        cols = zeros(res,res);
        cols(rows < 0) = cols_02(rows < 0);
        cols(rows > 0) = cols_01(rows > 0);
        
        im_03(:,:,n+1) = uint8(((cols / 2) .* window + 0.5) .* 255);
    end

    %% Sequences
    
    blank_im = 128 .* ones(res,res);
    
    im_01 = uint8(cat(3, blank_im, im_01));
    im_02 = uint8(cat(3, blank_im, im_02));
    im_03 = uint8(cat(3, blank_im, im_03));
    
    stim_frame  = 1/flicker_freq;
    
    if mod(frame_rate, flicker_freq)
        error('Stimulus frame is not a whole multiple of refresh rate')
    end
    
    total_sec = n_repeats * (period * 2) + period;
    n_frames = total_sec * flicker_freq;
    seqtiming = (0:n_frames-1) .* stim_frame;
    
    
    stimulus.cmap = gray(256);
    stimulus.srcRect = [0 0 res res];
    stimulus.destRect = [128 0 896 786];
    stimulus.seqtiming = seqtiming';
    
    to_add = mod(2:n_frames+(2-1),2)';
    
    for nn = 1:n_stim
        
        % fixation dot sequence
        fixSeq = ones(n_fix_changes,1)*round(rand(1,ceil(n_frames/n_fix_changes)));
        fixSeq = fixSeq(:)+1;
        stimulus.fixSeg = fixSeq(1:n_frames);
        
        
        eval(['stimulus.images = uint8(im_0' num2str(nn) ');']);
        blank_period = zeros(period * flicker_freq,1);
        
        if nn == 1
            tmp = repmat([1 0],1, n_repeats);
            tmp = tmp(ones(period*flicker_freq,1),:);
            tmp = [blank_period; tmp(:)];
            tmp(tmp~=0) = tmp(tmp~=0) + floor(rand(sum(tmp~=0),1) .* num_examples);
        else
            
            tmp = repmat([1 0],1, n_repeats) .* (1:n_repeats*2);
            tmp = tmp(ones(period*flicker_freq,1),:);
            tmp = [blank_period; tmp(:)];
            
            tmp(tmp~=0) = tmp(tmp~=0) + to_add(tmp~=0);
            
        end
        
        stimulus.seq = tmp+1;
        save(fullfile(save_path, sprintf('stimulus_fmri_carrier_%s_run%d',im_labels{nn},run)));
        
    end
end



