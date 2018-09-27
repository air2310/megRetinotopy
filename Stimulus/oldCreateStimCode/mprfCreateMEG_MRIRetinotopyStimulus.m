% MEG RETINOTOPY
%   MAKE CARRIER PATTERNS
%
% GOAL: Make multiple exemplars or several spatial patterns for testing
% their effectiveness in driving BOLD responses in visual cortex:
%
% (1) Narrow-band curvy patterns (ala SOC paper)
% (2) Dartboards
% (3) Checkerboards

% Assumes vistadisp is on you Matlab path and you're on the meg_compatible
% branch. Note to self: use git checkout master and git checkout
% meg_compatible
%
% addpath(genpath('~/Documents/MATLAB/toolboxes/vistadisp'));
%
% % Needs psychtoolbox 3
% tbUse('psychtoolbox-3');

%% In general: check wich has the smaller screen, relative to the subject's position:
cal_mri = 'CBI_NYU_projector';
cal_meg = 'meg_lcd';

w_pattern = 'checkerboard';
flicker_freqs = 10;

d_mri   = loadDisplayParams(cal_mri);
d_meg   = loadDisplayParams(cal_meg);
d_meg.dimensions = [23.9 16.6];
d_meg.distance = 45;

mri_res = min(d_mri.numPixels);
meg_res = min(d_meg.numPixels);

mri_fr  = d_mri.frameRate;
meg_fr  = d_meg.frameRate;

d_mri.pixelSize = min((d_mri.dimensions./d_mri.numPixels)); % cm per pixel
d_meg.pixelSize = min((d_meg.dimensions./d_meg.numPixels)); % cm per pixel

mri_scr_size = 2*pix2angle(d_mri,mri_res/2);
meg_scr_size = 2*pix2angle(d_meg,meg_res/2);

% This is the visual angle that we need to use for both experiments:
cur_screen_size = min([mri_scr_size meg_scr_size]);

mri_scale = (cur_screen_size / mri_scr_size);
meg_scale = (cur_screen_size / meg_scr_size);

%% Options:
bug_me_with_dialog = false;

if bug_me_with_dialog
    path = uigetfile('.mat','Select file to load stimulus parameters from');
    
else
    path = 0;
    
end

if path
    load(path)
    
    flicker_freqs = stimulus.opts.general.flicker_freqs;
    
    num_examples = stimulus.opts.patterns.num_examples;
    sparseness = stimulus.opts.patterns.sparseness;
    center_sf = stimulus.opts.patterns.center_sf;
    bandwidth = stimulus.opts.patterns.bandwidth;
    filter_size = stimulus.opts.patterns.filter_size;
    
    wegdes_dpc = stimulus.opts.dartboard.wedges_dpc;
    rings_dpc = stimulus.opts.dartboard.rings_dpc;
    n_wedges = stimulus.opts.dartboard.n_wedges;
    n_rings = stimulus.opts.dartboard.n_rings;
    
    rows_dpc = stimulus.opts.checkerboard.rows_dpc;
    cols_dpc = stimulus.opts.checkerboard.cols_dpc;
    n_rows = stimulus.opts.checkerboard.n_rows;
    n_cols = stimulus.opts.checkerboard.n_cols;
    
else
    fprintf('No file selected, loading defaults\n');
    
    if strcmpi(w_pattern, 'patterns')
        % Curvy patterns -> Check how this has to change to accommodate the
        % different screens:
        num_examples    = 10;  % Amount of patterns generated
        sparseness      = 1/10;% Sparseness of curvy patterns
        center_sf       = 3;    % cycles per degree
        bandwidth       = 2;    % Bandwidth of filter in octaves
        filter_size     = 31;   % size of filter in pixels
        
        % Addpaths
        addpath(genpath('~/matlab/git/soccode')) % Add paths and load the bpfilter and the base images
        addpath(genpath('~/matlab/git/knkutils'))
        addpath(genpath('~/matlab/git/matlabPyrTools'))
        
        w_filter_mri        = 'flt_dog_mri';
        w_filter_meg        = 'flt_dog_meg';

        
    elseif strcmpi(w_pattern,'dartboard');
        
        % Dartboard options
        wedges_dpc      = 2;    % i.e. degrees for one black/white repeat
        rings_dpc       = 2;
        
        
        n_wedges_meg        = meg_scr_size / wedges_dpc * 2; % Amount of wedges making up the dartboard
        n_rings_meg         = meg_scr_size / rings_dpc * 2;  % Amount of rings
        
        
        n_wedges_mri        = mri_scr_size / wedges_dpc * 2; % Amount of wedges making up the dartboard
        n_rings_mri         = mri_scr_size / rings_dpc * 2;  % Amount of rings
        
    elseif strcmpi(w_pattern,'checkerboard');
        % Checkerboard options
        rows_dpc        = 2;
        cols_dpc        = 2;
        
        n_rows_meg      = meg_scr_size / rows_dpc * 2; % Amount of checks along the board vertically
        n_cols_meg      = meg_scr_size / cols_dpc * 2; % Amount of checks along the board horizontally
        
        n_rows_mri      = mri_scr_size / rows_dpc * 2;
        n_cols_mri      = mri_scr_size / cols_dpc * 2;
        
    else
        error('Unknown pattern')
        
        
    end
    
end

%% Retinotopy options:
TR              = 1.5;
step_time_mri  = TR;
step_time_meg  = TR;

outer_rad       = cur_screen_size / 2;
flicker_freq    = 12;
n_fix_changes   = 20;   % amount of fixation color changes

stim_frame      = 1/flicker_freq;

mri_step_duration   = TR;
meg_step_duration   = TR;

bar_width           = cur_screen_size / 10;
step_size           = bar_width / 2;
add_coverage        = bar_width / 2;
low_x               = -(outer_rad+add_coverage) : step_size : outer_rad - add_coverage;
high_x              = low_x + bar_width;

orientations = (0:45:360)./360*(2*pi); % degrees -> rad
orientations = orientations([1 6 3 8 5 2 7 4]);


%% General stuff

[Xmeg, Ymeg] = meshgrid(-(meg_res-1)/2:(meg_res-1)/2);
[THmeg, Rmeg] = cart2pol(Xmeg,Ymeg);
THmeg(THmeg < 0) = THmeg(THmeg < 0) + 2 * pi;

[Xmri, Ymri] = meshgrid(-(mri_res-1)/2:(mri_res-1)/2);
[THmri, Rmri] = cart2pol(Xmri,Ymri);
THmri(THmri < 0) = THmri(THmri < 0) + 2 * pi;


Rmri_deg = pix2angle(d_mri,Rmri);
Rmeg_deg = pix2angle(d_meg,Rmeg);

Xmri_deg = pix2angle(d_mri,Xmri);
Xmeg_deg = pix2angle(d_meg,Xmeg);

Ymri_deg = pix2angle(d_mri,Ymri);
Ymeg_deg = pix2angle(d_meg,Ymeg);

% or
% Rmri_deg = Rmri .* pix2angle(d_mri,1);
% Rmeg_deg = Rmeg .* pix2angle(d_meg,1);
% Xmri_deg = Xmri .* pix2angle(d_mri,1);
% Xmeg_deg = Xmeg .* pix2angle(d_meg,1);
% Ymri_deg = Ymri .* pix2angle(d_mri,1);
% Ymeg_deg = Ymeg .* pix2angle(d_meg,1);

try
    [flt_cos_mri, flt_dog_mri] = mprfComputeAndStoreBPFilter('mri_filter');
    
catch
    [flt_cos_mri, flt_dog_mri] = mprfComputeAndStoreBPFilter([], center_sf, bandwidth, filter_size, mri_res * mri_scale, cur_screen_size, 'mri_filter',false);
    
end


try
    [flt_cos_meg, flt_dog_meg] = mprfComputeAndStoreBPFilter('meg_filter');
    
catch
    [flt_cos_meg, flt_dog_meg] = mprfComputeAndStoreBPFilter([], center_sf, bandwidth, filter_size, meg_res * meg_scale, cur_screen_size, 'meg_filter',false);
    
end

mri_pattern_images = nan(mri_res, mri_res, num_examples);
meg_pattern_images = nan(mri_res, mri_res, num_examples);

for n = 1:num_examples
    mri_pattern_images(:,:,n) = createPatternStimulus([mri_res, mri_res], sparseness, eval(w_filter_mri));
    meg_pattern_images(:,:,n) = createPatternStimulus([meg_res, meg_res], sparseness, eval(w_filter_meg));
    
end

mri_pattern_images = uint8(max(mri_pattern_images + 0.5,0) .*255);
meg_pattern_images = uint8(max(meg_pattern_images + 0.5,0) .*255);

im_mri_01 = uint8(zeros(mri_res, mri_res, flicker_freq * step_time_mri * length(orientations) * length(low_x)));
im_meg_01 = uint8(zeros(meg_res, meg_res, flicker_freq * step_time_meg * length(orientations) * length(low_x)));

im_mri_02 = uint8(zeros(mri_res, mri_res, 2 * length(orientations) * length(low_x)));
im_meg_02 = uint8(zeros(meg_res, meg_res, 2 * length(orientations) * length(low_x)));

im_mri_03 = uint8(zeros(mri_res, mri_res, 2 * length(orientations) * length(low_x)));
im_meg_03 = uint8(zeros(meg_res, meg_res, 2 * length(orientations) * length(low_x)));


cur_im = 0;
cur_im2 = 0;
cur_im3 = 0;

back_ground_mri = uint8(zeros(mri_res, mri_res) + 128);
back_ground_meg = uint8(zeros(meg_res, meg_res) + 128);



wedges_mri = sign(2*round((sin(THmri .* (n_wedges_mri/2))+1)/2)-1);
rings_01 = sign(2*round((sin( (Rmri./((mri_res-1)/2)).*2*pi .* (n_rings_mri/2))+1)/2)-1);
rings_02 = sign(2*round((sin( (Rmri./((mri_res-1)/2)).*2*pi .* (n_rings_mri/2) +pi)+1)/2)-1);

rings_mri_01 = zeros(mri_res,mri_res);
rings_mri_01(wedges_mri < 0) = rings_01(wedges_mri < 0);
rings_mri_01(wedges_mri > 0) = rings_02(wedges_mri > 0);

rings_mri_02 = zeros(mri_res,mri_res);
rings_mri_02(wedges_mri < 0) = rings_02(wedges_mri < 0);
rings_mri_02(wedges_mri > 0) = rings_01(wedges_mri > 0);



wedges_meg = sign(2*round((sin(THmeg .* (n_wedges_meg/2))+1)/2)-1);
rings_01 = sign(2*round((sin( (Rmeg./((meg_res-1)/2)).*2*pi .* (n_rings_meg/2))+1)/2)-1);
rings_02 = sign(2*round((sin( (Rmeg./((meg_res-1)/2)).*2*pi .* (n_rings_meg/2) +pi)+1)/2)-1);

rings_meg_01 = zeros(meg_res,meg_res);
rings_meg_01(wedges_meg < 0) = rings_01(wedges_meg < 0);
rings_meg_01(wedges_meg > 0) = rings_02(wedges_meg > 0);

rings_meg_02 = zeros(meg_res,meg_res);
rings_meg_02(wedges_meg < 0) = rings_02(wedges_meg < 0);
rings_meg_02(wedges_meg > 0) = rings_01(wedges_meg > 0);

first_idx = 1:2:size(im_meg_02,3);
second_idx = 2:2:size(im_meg_02,3);


fprintf('Creating %d images:',length(orientations) * length(low_x))

for n = 1:length(orientations)
    cur_or = orientations(n);
    
    newX_mri = Xmri .* cos(cur_or) - Ymri .* sin(cur_or);
    newX_meg = Xmeg .* cos(cur_or) - Ymeg .* sin(cur_or);
    
    newY_mri = Xmri .* sin(cur_or) + Ymri .* cos(cur_or);
    newY_meg = Xmeg .* sin(cur_or) + Ymeg .* cos(cur_or);
    
    newX_mri_deg = Xmri_deg .* cos(cur_or) - Ymri_deg .* sin(cur_or);
    newX_meg_deg = Xmeg_deg .* cos(cur_or) - Ymeg_deg .* sin(cur_or);
    
    
    rows_mri = sign(2*round((sin((newY_mri./((mri_res-1)/2)).*2*pi .* (n_rows_mri/4))+1)/2)-1);
    cols_mri_01 = sign(2*round((sin( (newX_mri./((mri_res-1)/2)).*2*pi .* (n_cols_mri/4))+1)/2)-1);
    cols_mri_02 = sign(2*round((sin( (newX_mri./((mri_res-1)/2)).*2*pi .* (n_cols_mri/4) +pi)+1)/2)-1);
    
    cols_01_mri = zeros(mri_res,mri_res);
    cols_01_mri(rows_mri < 0) = cols_mri_01(rows_mri < 0);
    cols_01_mri(rows_mri > 0) = cols_mri_02(rows_mri > 0);
    
    
    cols_02_mri = zeros(mri_res,mri_res);
    cols_02_mri(rows_mri < 0) = cols_mri_02(rows_mri < 0);
    cols_02_mri(rows_mri > 0) = cols_mri_01(rows_mri > 0);
    
    
    rows_meg = sign(2*round((sin((newY_meg./((meg_res-1)/2)).*2*pi .* (n_rows_meg/4))+1)/2)-1);
    cols_meg_01 = sign(2*round((sin( (newX_meg./((meg_res-1)/2)).*2*pi .* (n_cols_meg/4))+1)/2)-1);
    cols_meg_02 = sign(2*round((sin( (newX_meg./((meg_res-1)/2)).*2*pi .* (n_cols_meg/4) +pi)+1)/2)-1);
    
    cols_01_meg = zeros(meg_res,meg_res);
    cols_01_meg(rows_meg < 0) = cols_meg_01(rows_meg < 0);
    cols_01_meg(rows_meg > 0) = cols_meg_02(rows_meg > 0);
    
    
    cols_02_meg = zeros(meg_res,meg_res);
    cols_02_meg(rows_meg < 0) = cols_meg_02(rows_meg < 0);
    cols_02_meg(rows_meg > 0) = cols_meg_01(rows_meg > 0);
    
    
    for nn = 1:length(low_x)
        cur_im = cur_im + 1;
        
        if mod(cur_im, 10)
        else
            fprintf('.%d.',cur_im)
        end
        
        window_mri = ( (newX_mri_deg>=low_x(nn) & newX_mri_deg<=high_x(nn)) & Rmri_deg<outer_rad);
        window_meg = ( (newX_meg_deg>=low_x(nn) & newX_meg_deg<=high_x(nn)) & Rmeg_deg<outer_rad);
        %% Curvy patterns
        
        for nnn = 1:flicker_freq * step_time_mri
            idx = ceil(rand * num_examples);
            
            tmp = back_ground_mri;
            tmp2 = mri_pattern_images(:,:,idx);
            tmp(window_mri) = tmp2(window_mri);
            
            cur_im2 = cur_im2+1;
            im_mri_01(:,:,cur_im2) = tmp;
            
        end
        
        for nnn = 1:flicker_freq * step_time_meg
            cur_im3 = cur_im3+1;
            
            tmp = back_ground_meg;
            tmp2 = meg_pattern_images(:,:,idx);
            tmp(window_mri) = tmp2(window_mri);
            
            im_meg_01(:,:,cur_im3) = tmp;
            
        end
        
        %% Dartboard
        
        im_mri_02(:,:,first_idx(cur_im)) = uint8(((rings_mri_01 / 2) .* window_mri + 0.5) .* 255);
        im_mri_02(:,:,second_idx(cur_im)) = uint8(((rings_mri_02 / 2) .* window_mri + 0.5).*255);
        
        im_mri_03(:,:,first_idx(cur_im)) = uint8(((cols_01_mri / 2) .* window_mri + 0.5) .* 255);
        im_mri_03(:,:,second_idx(cur_im)) = uint8(((cols_02_mri / 2) .* window_mri + 0.5) .* 255);
        
        %% Checkerboard
        
        
        
        im_meg_02(:,:,first_idx(cur_im)) = uint8(((rings_meg_01 / 2) .* window_meg + 0.5) .* 255);
        im_meg_02(:,:,second_idx(cur_im)) = uint8(((rings_meg_01 / 2) .* window_meg + 0.5).*255);
        
        im_meg_03(:,:,first_idx(cur_im)) = uint8(((cols_01_meg / 2) .* window_meg + 0.5) .* 255);
        im_meg_03(:,:,second_idx(cur_im)) = uint8(((cols_01_meg / 2) .* window_meg + 0.5) .* 255);
        
        
        
        
        
    end
    
end

fprintf('.Done.\n');

return

figure;
for n = 1:size(im_mri_02,3)
    imshow(im_mri_03(:,:,n));
    drawnow;
    
    
end


%% Sequence for checks and dartboard

blank_im_mri = uint8(ones(mri_res, mri_res) .* 128);
blank_im_meg = uint8(ones(meg_res, meg_res) .* 128);


im_mri_01 = cat(3,blank_im_mri, im_mri_01);
im_meg_01 = cat(3,blank_im_meg, im_meg_01);

im_mri_02 = cat(3,blank_im_mri, im_mri_02);
im_meg_02 = cat(3,blank_im_meg, im_meg_02);

im_mri_03 = cat(3,blank_im_mri, im_mri_03);
im_meg_03 = cat(3,blank_im_meg, im_meg_03);


total_sec_mri = step_time_mri * length(orientations) * length(low_x);
total_sec_meg = step_time_meg * length(orientations) * length(low_x);

n_frames_mri = total_sec_mri * flicker_freq;
n_frames_meg = total_sec_meg * flicker_freq;

seqtiming_mri = (0:n_frames_mri-1) .* stim_frame;
seqtiming_meg = (0:n_frames_meg-1) .* stim_frame;

n_im_per_step_mri =  step_time_mri * flicker_freq /2;
n_im_per_step_meg =  step_time_meg * flicker_freq /2;

n_images_mri = size(im_mri_02,3)-1;
n_images_meg = size(im_meg_02,3)-1;


aa_mri = mod(2:n_frames_mri+(2-1),2)+1;
aa_meg = mod(2:n_frames_meg+(2-1),2)+1;

n_steps = total_sec_mri / step_time_mri

bb_mri = repmat(0:2:(2*n_steps)-1,(n_frames_mri / (n_images_mri/2)),1);
bb_meg = repmat(0:2:(2*n_steps)-1,(n_frames_meg / (n_images_meg/2)),1);

cc_mri = bb_mri(:);
cc_meg = bb_meg(:);

seq_mri = aa_mri(:) + cc_mri + 1;
seq_meg = aa_meg(:) + cc_meg + 1;

seqtiming_mri = [0:length(seq_mri)-1] .* stim_frame;
seqtiming_meg = [0:length(seq_meg)-1] .* stim_frame;
fixSeq_mri = ones(n_fix_changes,1)*round(rand(1,ceil(n_frames_mri/n_fix_changes)));
fixSeq_mri = fixSeq_mri(:)+1;

fixSeq_meg = ones(n_fix_changes,1)*round(rand(1,ceil(n_frames_meg/n_fix_changes)));
fixSeq_meg = fixSeq_meg(:)+1;



stimulus.seq = seq_mri;
stimulus.seqtiming = seqtiming_mri;
stimulus.images = im_mri_03;
stimulus.cmap = gray(256);
stimulus.fixSeq = fixSeq_mri(1:n_frames_mri);
stimulus.srcRect = [0 0 mri_res mri_res];
stimulus.destRect = [128 0 896 786];


save(fullfile('Stimuli','MRI_ret_stim_try'),'stimulus');
% stimulus.trigSeq
% stimulus.diodeSeg


figure;
colormap gray
for n = 1:length(dd_mri)
    imagesc(im_mri_02(:,:,dd_mri(n)));
    drawnow;
    pause(1/12)
    
    
end














































































