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
% Needs vistadisp (Winawerlab fork)
% addpath(genpath('~/matlab/git/forks/vistadisp'));
%
% Needs psychtoolbox 3
% tbUse('psychtoolbox-3'); % or
% addpath(genpath('/Applications/Psychtoolbox'))

%% In general: check wich has the smaller screen, relative to the subject's position:
cal_mri = 'CBI_Propixx';
cal_meg = 'meg_lcd';

mri_save_path = '/Volumes/server/Projects/MEG/Retinotopy/Stimuli/20180909/Retinotopy/MRI';
meg_save_path = '/Volumes/server/Projects/MEG/Retinotopy/Stimuli/20180909/Retinotopy/MEG';

stim_size = 'max'; % 'max' or the degrees visual angle it should be (number, not char)
nruns = 20;

w_pattern = 'checkerboard';
flicker_freq = 10;

d_mri   = loadDisplayParams(cal_mri);
d_meg   = loadDisplayParams(cal_meg);
% d_meg.dimensions = [23.9 16.6];
% d_meg.dimensions = 
% d_meg.distance = 45;

mri_res = min(d_mri.numPixels);
meg_res = min(d_meg.numPixels);

mri_fr  = d_mri.frameRate;
meg_fr  = d_meg.frameRate;

d_mri.pixelSize = min((d_mri.dimensions./d_mri.numPixels)); % cm per pixel
d_meg.pixelSize = min((d_meg.dimensions./d_meg.numPixels)); % cm per pixel

mri_scr_size = 2*pix2angle(d_mri,mri_res/2);
meg_scr_size = 2*pix2angle(d_meg,meg_res/2);

if ischar('max')
    if strcmpi(stim_size,'max')
        
        % This is the visual angle that we need to use for both experiments if fullscreen:
        cur_screen_size = min([mri_scr_size meg_scr_size]);
        
    else
        error('Unknown stimulus size')
        
    end
elseif isnumeric(stim_size)
    % We can set it manually as well:
    cur_screen_size = stim_size;
    tmp = [mri_scr_size meg_scr_size] .* (cur_screen_size / max([mri_scr_size meg_scr_size]));
    mri_scr_size = tmp(1);
    meg_scr_size = tmp(2);
    
end

mri_scale = (cur_screen_size / mri_scr_size);
meg_scale = (cur_screen_size / meg_scr_size);

%% Options:


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
    
    
elseif strcmpi(w_pattern,'dartboard')
    
    % Dartboard options
    n_wedges        = 20; % Amount of wedges making up the dartboard
    rings_dpc       = 2;  % i.e. degrees visual angle for one black/white repeat
    
    
    n_wedges_meg        = n_wedges;  % Amount of wedges making up the dartboard
    n_rings_meg         = meg_scr_size / rings_dpc * 2;  % Amount of rings
    
    
    n_wedges_mri        = n_wedges; % Amount of wedges making up the dartboard
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


%% Retinotopy options:
TR              = 1.5;
step_time_mri  = TR;
step_time_meg  = TR;

blank_time_mri = 15 * TR;
blink_time_meg = 2 * step_time_meg;
blank_time_meg = 3 * step_time_meg; % Also time for blinking? Better not, as we really want to use these as a baseline measurement.
% Preferably, blink period directly after stimulus presentation that gets
% descarded, than a blank period

n_repeats_first_stim_time_meg = 1; % As it appears to be common practice to discard the first second/epoch of a new stimulus sequence
% This option will allow you to repeat the first stimulus position of a
% bar sweep.

blink_sq_size   = 2;    % Degrees of blink square (width)
blink_size      = round(2*angle2pix(d_meg,blink_sq_size/2));

outer_rad       = cur_screen_size / 2;
n_fix_changes   = 20;   % amount of fixation color changes

stim_frame      = 1/flicker_freq;

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

for run = 1:nruns
    
    if strcmpi(w_pattern, 'patterns')
        
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
        
        mri_raw_images = nan(mri_res, mri_res, num_examples);
        meg_raw_images = nan(mri_res, mri_res, num_examples);
        
        for n = 1:num_examples
            mri_raw_images(:,:,n) = createPatternStimulus([mri_res, mri_res], sparseness, eval(w_filter_mri));
            meg_raw_images(:,:,n) = createPatternStimulus([meg_res, meg_res], sparseness, eval(w_filter_meg));
            
        end
        
        mri_raw_images = uint8(max(mri_raw_images + 0.5,0) .*255);
        meg_raw_images = uint8(max(meg_raw_images + 0.5,0) .*255);
        
        
        im_mri = uint8(zeros(mri_res, mri_res, flicker_freq * step_time_mri * length(orientations) * length(low_x)));
        im_meg = uint8(zeros(meg_res, meg_res, flicker_freq * step_time_meg * length(orientations) * length(low_x)));
        
        
    elseif strcmpi(w_pattern, 'dartboard')
        
        
        
        
        wedges_mri = sign(2*round((sin(THmri .* (n_wedges_mri/2))+1)/2)-1);
        rings_01 = sign(2*round((sin( (Rmri./((mri_res-1)/2)).*2*pi .* (n_rings_mri/2))+1)/2)-1);
        rings_02 = sign(2*round((sin( (Rmri./((mri_res-1)/2)).*2*pi .* (n_rings_mri/2) +pi)+1)/2)-1);
        
        mri_raw_images_01 = zeros(mri_res,mri_res);
        mri_raw_images_01(wedges_mri < 0) = rings_01(wedges_mri < 0);
        mri_raw_images_01(wedges_mri > 0) = rings_02(wedges_mri > 0);
        
        mri_raw_images_02 = zeros(mri_res,mri_res);
        mri_raw_images_02(wedges_mri < 0) = rings_02(wedges_mri < 0);
        mri_raw_images_02(wedges_mri > 0) = rings_01(wedges_mri > 0);
        
        
        
        wedges_meg = sign(2*round((sin(THmeg .* (n_wedges_meg/2))+1)/2)-1);
        rings_01 = sign(2*round((sin( (Rmeg./((meg_res-1)/2)).*2*pi .* (n_rings_meg/2))+1)/2)-1);
        rings_02 = sign(2*round((sin( (Rmeg./((meg_res-1)/2)).*2*pi .* (n_rings_meg/2) +pi)+1)/2)-1);
        
        meg_raw_images_01 = zeros(meg_res,meg_res);
        meg_raw_images_01(wedges_meg < 0) = rings_01(wedges_meg < 0);
        meg_raw_images_01(wedges_meg > 0) = rings_02(wedges_meg > 0);
        
        meg_raw_images_02 = zeros(meg_res,meg_res);
        meg_raw_images_02(wedges_meg < 0) = rings_02(wedges_meg < 0);
        meg_raw_images_02(wedges_meg > 0) = rings_01(wedges_meg > 0);
        
        
        
        im_mri = uint8(zeros(mri_res, mri_res, 2 * length(orientations) * length(low_x)));
        im_meg = uint8(zeros(meg_res, meg_res, 2 * length(orientations) * length(low_x)));
        
        
    elseif strcmpi(w_pattern,'checkerboard');
        
        
        im_mri = uint8(zeros(mri_res, mri_res, 2 * length(orientations) * length(low_x)));
        im_meg = uint8(zeros(meg_res, meg_res, 2 * length(orientations) * length(low_x)));
        
        
    else
        error('Unknown pattern')
        
    end
    
    
    
    cur_im = 0;
    cur_im2 = 0;
    cur_im3 = 0;
    
    back_ground_mri = uint8(zeros(mri_res, mri_res) + 128);
    back_ground_meg = uint8(zeros(meg_res, meg_res) + 128);
    
    all_orientations = nan(1,size(im_mri,3));
    
    
    
    first_idx = 1:2:size(im_mri,3);
    second_idx = 2:2:size(im_mri,3);
    
    
    fprintf('Creating %d images:',length(orientations) * length(low_x))
    
    for n = 1:length(orientations)
        cur_or = orientations(n);
        newX_mri_deg = Xmri_deg .* cos(cur_or) - Ymri_deg .* sin(cur_or);
        newX_meg_deg = Xmeg_deg .* cos(cur_or) - Ymeg_deg .* sin(cur_or);
        rph1 = rand .* (2*pi);
        rph2 = rand .* (2*pi);
        rph3 = rand .* (2*pi);
        rph4 = rand .* (2*pi);
        
        
        if strcmpi(w_pattern, 'checkerboard')
            newX_mri = Xmri .* cos(cur_or) - Ymri .* sin(cur_or);
            newX_meg = Xmeg .* cos(cur_or) - Ymeg .* sin(cur_or);
            
            newY_mri = Xmri .* sin(cur_or) + Ymri .* cos(cur_or);
            newY_meg = Xmeg .* sin(cur_or) + Ymeg .* cos(cur_or);
            
            
            
            rows_mri = sign(2*round((sin((newY_mri./((mri_res-1)/2)).*2*pi .* (n_rows_mri/4) + rph1)+1)/2)-1);
            cols_mri_01 = sign(2*round((sin( (newX_mri./((mri_res-1)/2)).*2*pi .* (n_cols_mri/4) + rph2)+1)/2)-1);
            cols_mri_02 = sign(2*round((sin( (newX_mri./((mri_res-1)/2)).*2*pi .* (n_cols_mri/4) +pi + rph2)+1)/2)-1);
            
            mri_raw_images_01 = zeros(mri_res,mri_res);
            mri_raw_images_01(rows_mri < 0) = cols_mri_01(rows_mri < 0);
            mri_raw_images_01(rows_mri > 0) = cols_mri_02(rows_mri > 0);
            
            
            mri_raw_images_02 = zeros(mri_res,mri_res);
            mri_raw_images_02(rows_mri < 0) = cols_mri_02(rows_mri < 0);
            mri_raw_images_02(rows_mri > 0) = cols_mri_01(rows_mri > 0);
            
            
            rows_meg = sign(2*round((sin((newY_meg./((meg_res-1)/2)).*2*pi .* (n_rows_meg/4) + rph3)+1)/2)-1);
            cols_meg_01 = sign(2*round((sin( (newX_meg./((meg_res-1)/2)).*2*pi .* (n_cols_meg/4) + rph4)+1)/2)-1);
            cols_meg_02 = sign(2*round((sin( (newX_meg./((meg_res-1)/2)).*2*pi .* (n_cols_meg/4) +pi + rph4)+1)/2)-1);
            
            meg_raw_images_01 = zeros(meg_res,meg_res);
            meg_raw_images_01(rows_meg < 0) = cols_meg_01(rows_meg < 0);
            meg_raw_images_01(rows_meg > 0) = cols_meg_02(rows_meg > 0);
            
            
            meg_raw_images_02 = zeros(meg_res,meg_res);
            meg_raw_images_02(rows_meg < 0) = cols_meg_02(rows_meg < 0);
            meg_raw_images_02(rows_meg > 0) = cols_meg_01(rows_meg > 0);
        end
        
        for nn = 1:length(low_x)
            cur_im = cur_im + 1;
            all_orientations(first_idx(cur_im)) = cur_or / pi * 4; % multiples of 45
            all_orientations(second_idx(cur_im)) = cur_or / pi * 4; % multiples of 45
            
            if mod(cur_im, 10)
            else
                fprintf('.%d.',cur_im)
            end
            
            window_mri = ( (newX_mri_deg>=low_x(nn) & newX_mri_deg<=high_x(nn)) & Rmri_deg<outer_rad);
            window_meg = ( (newX_meg_deg>=low_x(nn) & newX_meg_deg<=high_x(nn)) & Rmeg_deg<outer_rad);
            
            if strcmpi(w_pattern,'patterns')
                %% Curvy patterns
                
                for nnn = 1:flicker_freq * step_time_mri
                    idx = ceil(rand * num_examples);
                    
                    tmp = back_ground_mri;
                    tmp2 = mri_raw_images(:,:,idx);
                    tmp(window_mri) = tmp2(window_mri);
                    
                    cur_im2 = cur_im2+1;
                    im_mri(:,:,cur_im2) = tmp;
                    
                end
                
                for nnn = 1:flicker_freq * step_time_meg
                    cur_im3 = cur_im3+1;
                    
                    tmp = back_ground_meg;
                    tmp2 = meg_raw_images(:,:,idx);
                    tmp(window_meg) = tmp2(window_meg);
                    
                    im_meg(:,:,cur_im3) = tmp;
                    
                end
                
            elseif strcmpi(w_pattern,'checkerboard') || ...
                    strcmpi(w_pattern,'dartboard')
                
                im_mri(:,:,first_idx(cur_im)) = uint8(((mri_raw_images_01 / 2) .* window_mri + 0.5) .* 255);
                im_mri(:,:,second_idx(cur_im)) = uint8(((mri_raw_images_02 / 2) .* window_mri + 0.5).*255);
                
                im_meg(:,:,first_idx(cur_im)) = uint8(((meg_raw_images_01 / 2) .* window_meg + 0.5) .* 255);
                im_meg(:,:,second_idx(cur_im)) = uint8(((meg_raw_images_02 / 2) .* window_meg + 0.5) .* 255);
                
            else
                error('Unknown pattern')
                
            end
            
            
            
        end
        
    end
    
    fprintf('.Done.\n');
    
    %
    % figure;
    % for n = 1:size(im_mri,3)
    %     imshow(im_meg(:,:,n));
    %     drawnow;
    %
    %
    % end
    
    %% Sequence for checks and dartboard
    
    % Create blank image (gray / mean luminance)
    blank_im_mri = uint8(ones(mri_res, mri_res) .* 128);
    blank_im_meg = uint8(ones(meg_res, meg_res) .* 128);
    
    mri_ori_idx = [1 3 5 7];
    meg_ori_idx = 1:7;
    
    % When to add blank periods in sequence?
    blanks_after_these_orientations_mri = orientations(mri_ori_idx) ./ pi * 4;
    blanks_after_these_orientations_meg = orientations(meg_ori_idx) ./ pi * 4;
    
    % Create blink image (a mean luminance screen with a black square in
    % the middle)
    blink_im = zeros(1,meg_res);
    blink_im((meg_res/2 - round(blink_size/2)) : (meg_res/2 - round(blink_size/2)) + blink_size-1) = 1;
    blink_im = -1 .* (blink_im' * blink_im);
    blink_im = blink_im .*128 + 128;
    
	% Concatenate all images, both blank and bar (and blink block for MEG)
    im_mri_all = cat(3,blank_im_mri, im_mri);
    im_meg_all = cat(3,blank_im_meg, blink_im, im_meg);
    
    % Number of total frames in seconds: MRI: 342 (TR * number of orientations *
    % step size --dependent on the height of the smaller screen size + TR*blank time)
    total_sec_mri = step_time_mri * length(orientations) * length(low_x) + (length(mri_ori_idx)* blank_time_mri);
    total_sec_meg = step_time_meg * (length(orientations) * (length(low_x)+ n_repeats_first_stim_time_meg)) + (length(meg_ori_idx)* (blink_time_meg + blank_time_meg));
    
    % Number of stimulus frames in seconds: MRI: 252 (TR * number of orientations *
    % step size (dependent on the height of the smaller screen size)
    stim_sec_mri = step_time_mri * length(orientations) * length(low_x);
    stim_sec_meg = step_time_meg * (length(orientations) * (length(low_x) + n_repeats_first_stim_time_meg));
    
    % Number of total frames: MRI: 3420 (taking into account the 10 Hz flicker rate)  
    n_frames_mri = total_sec_mri * flicker_freq;
    n_frames_meg = total_sec_meg * flicker_freq;
    
    % Number of stimulus frames: MRI: 2520 (taking into account the 10 Hz flicker rate)  
    n_stim_frames_mri = stim_sec_mri * flicker_freq;
    n_stim_frames_meg = stim_sec_meg * flicker_freq;
    
    % Number of blank frames in between the bar sweeps: MRI: 255 (taking into account the 10 Hz flicker rate)  
    blank_frames_pp_mri = blank_time_mri * flicker_freq;
    blank_frames_pp_meg = blank_time_meg * flicker_freq;
    blink_frames_pp_meg = blink_time_meg * flicker_freq;
    
    n_images_mri = size(im_mri_all,3)-1;
    n_images_meg = size(im_meg_all,3)-2;
    %
    % n_im_per_step_mri =  step_time_mri * flicker_freq /2;
    % n_im_per_step_meg =  step_time_meg * flicker_freq /2;
    
    %% Compile MRI sequence
    
    % Use all frames to get an alternating sequence (1, 2) to get the flickering 
    aa_mri = mod(2:n_stim_frames_mri+(2-1),2)+1;
    
    % Get number of steps (frames) to get through all stimuli with this TR
    n_steps = stim_sec_mri / step_time_mri;
    
    % Put a stimulus frames into a sequence that shows a frame every other TR.
    % This sequence repeats 15 times since we have 15 conditions (stimulus sweeps and blanks)
    bb_mri = repmat(0:2:(2*n_steps)-1,(n_stim_frames_mri / (n_images_mri/2)),1);
    
    % Concatenate the alternating sequence of the 15 stimulus conditions
    cc_mri = bb_mri(:);
    
    % Merge the alternating (1,2) sequence with the stimulus condition
    % sequence
    im_seq_mri = aa_mri(:) + cc_mri+1;
    
    % Get nr of frames per pass (or sweep) of an oriented bar  (315 = 2520 / 8)
    mri_frames_per_pass = length(im_seq_mri) / length(orientations);
    
    % Get frame where each of the four blank periods will start
    blank_start_idx_mri = mri_ori_idx * mri_frames_per_pass + ...
        [1 blank_frames_pp_mri .* (1:length(mri_ori_idx)-1) + 1];
    
    % Get frame where each of the four blank periods will end
    blank_end_idx_mri = blank_start_idx_mri + blank_frames_pp_mri-1;
    
    % Get frame where each of the 5 corresponding stim periods to those blanks will start
    stim_start_idx = [0 blank_end_idx_mri] + 1;
    
    % Get frame where each of the 5 corresponding stim periods to those
    % blanks will end
    stim_end_idx = [blank_start_idx_mri - 1 n_frames_mri];
    
    % Preallocate space for the full mri sequence
    seq_mri = nan(1,n_frames_mri);
    
    % How many bar passes before a blank period?
    n_bar_passes = [mri_ori_idx(1) diff(mri_ori_idx) (length(orientations) - mri_ori_idx(end))];
    
    % Loop over all stimulus sweeps
    % Sequence is now: Stim, Blank, Stim, Stim, Blank, Stim, Stim, Blank,
    % Stim Stim, Blank, Stim. (One stim is 21 TRs, one Blank is 15 TRs) 
    start_frame = 1;
    for n = 1:length(stim_start_idx)
        s_idx = stim_start_idx(n);
        e_idx = stim_end_idx(n);
        n_frames = n_bar_passes(n) * mri_frames_per_pass;
        
        seq_mri(s_idx:e_idx) = im_seq_mri(start_frame:start_frame+n_frames-1);
        
        start_frame = start_frame + n_frames;
        
    end
    
    if sum(isnan(seq_mri)) == 4 * blank_frames_pp_mri;
        
    else
        error('Something went wrong')
    end
    
    % Add addition pre and post blank to mri sequence (8 trs)
    seq_mri(isnan(seq_mri)) = 1;
    n_tr_blank = 8;
    seq_mri = [ones(1, flicker_freq * TR * n_tr_blank) seq_mri ones(1, flicker_freq * TR * n_tr_blank)];
    
    seqtiming_mri = [0:length(seq_mri)-1] .* stim_frame;
    fixSeq_mri = ones(n_fix_changes,1)*round(rand(1,ceil(n_frames_mri/n_fix_changes)));
    fixSeq_mri = fixSeq_mri(:)+1;
    
    
    stimulus.seq = seq_mri;
    stimulus.seqtiming = seqtiming_mri;
    stimulus.images = im_mri_all;
    stimulus.cmap = gray(256);
    stimulus.fixSeq = fixSeq_mri(1:n_frames_mri);
    stimulus.srcRect = [0 0 mri_res mri_res];
    stimulus.destRect = [128 0 896 768];
    
    
%     save(fullfile(mri_save_path,sprintf('MRI_retinotopy_stimulus_run_%d',run)),'stimulus');
    clear stimulus
    
    
    % Compile MEG sequence
    aa_meg = mod(2:n_stim_frames_meg+(2-1),2)+1;
    n_steps = stim_sec_meg / step_time_meg;
    if n_repeats_first_stim_time_meg
        
        
        tmp = [zeros(1,n_repeats_first_stim_time_meg) 0:2:(2*(n_steps/length(orientations) - n_repeats_first_stim_time_meg))-1];
        tmp2 = repmat(tmp,1,length(orientations));
        tmp3 = repmat(0:length(low_x) * 2: (size(im_meg,3))-1, n_steps/length(orientations),1);
        bb_meg = repmat((tmp2(:) + tmp3(:))',flicker_freq * step_time_meg,1);
        % Need to repeat this with an extra column for every orientation:
        
    else
        bb_meg = repmat(0:2:(2*n_steps)-1,(n_stim_frames_meg / (n_images_meg/2)),1);
        
        
    end
    
    cc_meg = bb_meg(:);
    im_seq_meg = aa_meg(:) + cc_meg + 2;
    
    meg_frames_per_pass = length(im_seq_meg) / length(orientations);
    
    blank_start_idx_meg = meg_ori_idx * meg_frames_per_pass+ ...
        [1 (blank_frames_pp_meg + blink_frames_pp_meg) .* (1:length(meg_ori_idx)-1) + 1];
    
    blank_end_idx_meg = blank_start_idx_meg + blank_frames_pp_meg + blink_frames_pp_meg -1;
    
    stim_start_idx = [0 blank_end_idx_meg] + 1;
    stim_end_idx = [blank_start_idx_meg - 1 n_frames_meg];
    
    seq_meg = nan(1,n_frames_meg);
    n_bar_passes = [meg_ori_idx(1) diff(meg_ori_idx) (length(orientations) - meg_ori_idx(end))];
    
    start_frame = 1;
    for n = 1:length(stim_start_idx)
        s_idx = stim_start_idx(n);
        e_idx = stim_end_idx(n);
        n_frames = n_bar_passes(n) * meg_frames_per_pass;
        
        seq_meg(s_idx:e_idx) = im_seq_meg(start_frame:start_frame+n_frames-1);
        
        start_frame = start_frame + n_frames;
        
    end
    
    if sum(isnan(seq_meg)) == length(meg_ori_idx) * (blank_frames_pp_meg + blink_frames_pp_meg);
        
    else
        error('Something went wrong')
    end
    
    seq_meg(isnan(seq_meg)) = 1;
    all_triggers = [10 20 all_orientations+1]; % 10 = blink, 20 = blank orientations - 1 * 45 = orientation in degrees
    
    n_blink_frames = blink_time_meg ./ step_time_meg;
    n_blank_frames = blank_time_meg ./ step_time_meg;
    
    blink_idx = [zeros(1,length(low_x) + n_repeats_first_stim_time_meg) ones(1,n_blink_frames) zeros(1,n_blank_frames)];
    blink_idx_all = repmat(blink_idx,flicker_freq * step_time_meg,length(orientations));
    blink_idx_all = blink_idx_all(:);
    blink_idx_all = blink_idx_all(1:n_frames_meg);
    seq_meg(logical(blink_idx_all)) = 2;
    
    seqtiming_meg = [0:length(seq_meg)-1] .* stim_frame;
    trigger_idx = mod(seqtiming_meg, step_time_meg) == 0; % Send a trigger every step_time_interval
    trigger_seq_all = all_triggers(seq_meg);
    trigger_seq = trigger_seq_all .* trigger_idx;
    diode_seq = (trigger_seq ~=0);
    
    fixSeq_meg = ones(n_fix_changes,1)*round(rand(1,ceil(n_frames_meg/n_fix_changes)));
    fixSeq_meg = fixSeq_meg(:)+1;
    
    
    
    stimulus.seq = seq_meg;
    stimulus.seqtiming = seqtiming_meg;
    stimulus.images = im_meg_all;
    stimulus.cmap = gray(256);
    stimulus.fixSeq = fixSeq_meg(1:n_frames_meg);
    stimulus.srcRect = [0 0 meg_res meg_res];
    stimulus.destRect = [128 0 896 786];
    stimulus.trigSeq = trigger_seq;
    stimulus.diodeSeq = diode_seq;
    
    save(fullfile(meg_save_path,sprintf('MEG_retinotopy_stimulus_run_%d',run)),'stimulus');
    clear stimulus
end

return
figure;
colormap gray
for n = 1:length(seq_meg)
    imagesc(im_meg_all(:,:,seq_meg(n)));
    drawnow;
    pause(1/10)
    
end


%% To remember:
% adjust values on all_orientations, one orientation is zero now, i.e. no
% trigger?

% Use n_repeat_first_stim_time_meg to repeat first stim epoch for a certain
% fixed amount of times




























