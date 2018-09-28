% MEG RETINOTOPY
%   MAKE CARRIER PATTERNS FOR DIFFERENT FLICKER FREQUENCIES
%
% GOAL: Make multiple exemplars or several spatial patterns for testing
% their effectiveness in driving BOLD responses in visual cortex:
%
% (1) Narrow-band curvy patterns (ala SOC paper)
% (2) Dartboards
% (3) Checkerboards
%
% stimulus on/off periods can be varied independently. Make sure the
% spatial parameters of the eventual stimulus are the same as the one used
% in the spatial pattern experiment

% To do:
% make trigger sequence
% make diode sequence
% determine when subject can blink during blanks
% Look at:
% ExampleStimulusMEG.mat
% But also at:
%
%% Stimulus parameters:

% General
n_stim          = 3;    % one for every pattern
im_labels       = {'patterns','dartboard','checkerboard'}; % Names of patterns for saving
n_repeats       = 1;    % cycles per stimulus
n_fix_changes   = 20;   % amount of fixation color changes
n_runs          = 1;

on_period       = 12; % Seconds stim on/off
off_period      = 5; % Seconds stimulus off
blink_period    = 2; % of which this amount blink period
no_blink_period_01 = 1.5;
no_blink_period_02 = 1.5;
blink_sq_size   = 2;    % Degrees of blink square (width)

flicker_freqs   = [30 20 15 12 10 6 5 4 ]; %60 ./ (2 : 12); % Hz, i.e. images per second
save_path       = '/Volumes/server/Projects/MEG/Retinotopy/Stimuli/Flicker_separate_ff';

for n = 1:length(flicker_freqs)
    if any(mod([on_period off_period blink_period no_blink_period_01+no_blink_period_02],1/(flicker_freqs(n))))
        error('Periods are not consistent with flicker frequency %0.2f',flicker_freqs(n));
    end
    
end
% Screen options
cal = 'meg_lcd';
d = loadDisplayParams(cal);

% This is what we measured on 10 - 03 - 2017:
d_meg.dimensions = [23.9 16.6];
d_meg.distance  = 45;
d.pixelSize     = min(d.dimensions ./ d.numPixels);
res             = min(d.numPixels); % Size of stimulus in pixels. Stimulus will seen through circular window with this diameter
frame_rate      = d.frameRate;
screen_size     = 2*pix2angle(d,res/2);
blink_size      = round(2*angle2pix(d,blink_sq_size/2));

% Curvy patterns
num_examples    = 32;  % Amount of patterns generated
sparseness      = 1/10;% Sparseness of curvy patterns
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
opts.general.flicker_freqs = flicker_freqs;
opts.general.screen_size = screen_size;


% Store curvy pattern options:
opts.patterns.num_examples = num_examples;
opts.patterns.sparseness = sparseness;
opts.patterns.filter_center = center_sf;
opts.patterns.filter_bw = bandwidth;
opts.patterns.filter_size = filter_size;
opts.general.w_filter = w_filter;


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
    for cur_ff = 1:length(flicker_freqs);
        flicker_freq = flicker_freqs(cur_ff);
        
        
        try
            [flt_cos, flt_dog] = mprfComputeAndStoreBPFilter('pattern_filters');
            
        catch
            [flt_cos, flt_dog] = mprfComputeAndStoreBPFilter([], center_sf, bandwidth, filter_size, res, screen_size, 'pattern_filters',true);
            
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
        
        blink_im = zeros(1,res);
        blink_im((res/2 - round(blink_size/2)) : (res/2 - round(blink_size/2)) + blink_size-1) = 1;
        blink_im = -1 .* (blink_im' * blink_im);
        blink_im = blink_im .*128 + 128;
        
        im_01 = uint8(cat(3, blank_im, im_01, blink_im));
        im_02 = uint8(cat(3, blank_im, im_02, blink_im));
        im_03 = uint8(cat(3, blank_im, im_03, blink_im));
        
        stim_frame  = 1/flicker_freq;
        
        if mod(frame_rate, flicker_freq)
            error('Stimulus frame is not a whole multiple of refresh rate')
        end
        
        total_sec = n_repeats * (off_period + on_period);
        n_frames = total_sec * flicker_freq;
        seqtiming = (0:n_frames-1) .* stim_frame;
        
        
        stimulus.cmap = gray(256);
        stimulus.srcRect = [0 0 res res];
        stimulus.destRect = [128 0 896 786];
        stimulus.seqtiming = seqtiming';
        
        to_add = mod(2:n_frames+(2-1),2)';
        
        
        % First stim on frame, first stim off frame
        stim_on_off_frame = [1; zeros((on_period * flicker_freq) -1,1);1; zeros((off_period * flicker_freq)-1,1)];
        stim_on_off_frame = repmat(stim_on_off_frame ,1, n_repeats);
        stim_on_off_frame = stim_on_off_frame(:);
        
        % First stim on frame, first stim off frame
        stim_or_blank = [1; zeros((on_period * flicker_freq) -1,1);1; zeros((off_period * flicker_freq)-1,1)];
        stim_or_blank = repmat(stim_or_blank ,1, n_repeats);
        stim_or_blank = stim_or_blank(:);
        
        % Blink frames
%         blink_seq = [zeros((on_period + no_blink_period_01) * flicker_freq,1); ones(blink_period * flicker_freq,1);zeros(flicker_freq * no_blink_period_02,1)];
        blink_seq = [zeros((on_period) * flicker_freq,1); ones(blink_period * flicker_freq,1);zeros(flicker_freq * (no_blink_period_01 + no_blink_period_02),1)];

        blink_seq = repmat(blink_seq ,1, n_repeats+1);
        blink_seq = blink_seq(:);
        
%         trigSeq = [1; zeros((on_period * flicker_freq) -1,1);2; zeros(((no_blink_period_01+blink_period) * flicker_freq)-1,1); 2; zeros((no_blink_period_02 * flicker_freq)-1,1)];
         trigSeq = [flicker_freq; zeros(((on_period + blink_period) * flicker_freq) -1,1);2; zeros(((no_blink_period_01 + no_blink_period_02) * flicker_freq)-1,1)];

        trigSeq = repmat(trigSeq ,1, n_repeats);
        trigSeq = trigSeq(:);
        
        tmp = diff([0;blink_seq]) > 0;
        
        trigSeq(tmp(1:length(trigSeq))) = 3;
        
        
        
        for nn = 1:n_stim
            
            % fixation dot sequence
            fixSeq = ones(n_fix_changes,1)*round(rand(1,ceil(n_frames/n_fix_changes)));
            fixSeq = fixSeq(:)+1;
            stimulus.fixSeq = fixSeq(1:n_frames);
            
            
            eval(['stimulus.images = uint8(im_0' num2str(nn) ');']);
            num_im = size(stimulus.images,3);
            
            
            if nn == 1
                tmp = [ones(on_period * flicker_freq,1); zeros(off_period * flicker_freq,1)];
                tmp = repmat(tmp ,1, n_repeats);
                
                tmp = tmp(:);
                tmp(tmp~=0) = tmp(tmp~=0) + floor(rand(sum(tmp~=0),1) .* num_examples);
                
            else
                tmp = [ones(on_period * flicker_freq,1); zeros(off_period * flicker_freq,1)];
                tmp = bsxfun(@times, repmat(tmp ,1, n_repeats), (1:2:n_repeats*2));
                tmp = tmp(:);
                
                tmp(tmp~=0) = tmp(tmp~=0) + to_add(tmp~=0);
                
            end
            
            stimulus.trigSeq = trigSeq;
            stimulus.diodeSeq = stim_on_off_frame(:);
            stimulus.seq = tmp+1;
            stimulus.seq(logical(blink_seq(1:length(stimulus.seq)))) = num_im;
            
            stimulus.opts = opts;
            save(fullfile(save_path, sprintf('stimulus_fmri_carrier_%s_run%d_ff%d',im_labels{nn},run,flicker_freq)), 'stimulus');
            
        end
    end
end


return
% Check stuff:
figure;
colormap gray
for nn = 1:size(stimulus.seqtiming,1)
    imagesc(stimulus.images(:,:,stimulus.seq(nn)));
    drawnow;
    pause(diff(stimulus.seqtiming([1 2])));
    
    
end


% cond = 1;
%
% figure;
% colormap gray
% for nn = 1:size(all_cond,1)
%     imagesc(images(:,:,all_cond(nn,cond)));
%     drawnow;
%     pause(stim_frame);
%
%
% end



% for n = 1:2:N_blocks*2
%    tmp = repmat([n 0],1, n_repeats);
%    tmp = tmp(ones(period*flicker_freq,1),:);
%    tmp = tmp(:);
%    if n == 5
%        tmp(tmp~=0) = tmp(tmp~=0) + floor(rand(sum(tmp~=0),1) .* num_examples);
%
%    else
%        tmp(tmp~=0) = tmp(tmp~=0) + to_add(tmp~=0);
%
%    end
%    all_cond = [all_cond tmp];
%
% end






