%% s_BensonAtlasPRFPredictions
%
% This script is to make predicted responses for MEG sensors, given a 
% retinotopy bar stimulus, using the pRF of every vertex from the Benson et
% al. (2014) atlas of V1-V3. 
%
% This script relies on multiple preprocessing steps being complete:
% (1) subject has an autosegmented FS subject dir
% (2) subject has Benson atlas in the surf directory
% (3) subject has 'normal' FS surface (not downsampled) mesh in Brainstorm
% and Gain matrix computed from this FS surface
% (headmodel_surf_os_meg_02.mat)
% (4) MEG stimulus

% Script overview:
% 0. Set up: find paths and files, define parameters
% 1. Get Benson V1-V3 template that matches BS vertices (high res == FS
% vertices, low res == ~15,000 vertices)
% 2. Load templates
% 3. Get x,y coords from polar angle and eccen
% 4. Get RF from x, y, sigma and MEG stim aperture window
% 5. Compute predicted BOLD response for every vertex in V1-V3, where pRF has unit volume
% 6. Load Gain matrix and set dipole perpendicular to surface
% 7. Compute predicted MEG sensor responses and plot against data

%% 0. Set up

% Check if freesurfer matlab toolbox is added
if ~exist('MRIread') 
    addpath(genpath('/Applications/freesurfer/matlab'));
    addpath(genpath('/Applications/freesurfer/fsfast/toolbox'))
end

% Define params
subject     = 'wlsubj004';
bsProtocol  = 'MEG_Retinotopy';
highres     = true; % use 'normal' FS surface instead of downsampled surface

% Define directories
fsdir    = '/Volumes/server/Freesurfer_subjects'; % Freesurfer subjects 
bsDB     = '/Volumes/server/Projects/MEG/brainstorm_db'; % Brainstorm Data base
megRetDir   = '/Volumes/server/Projects/MEG/Retinotopy/'; % MEG retinotopy project on server



d = dir(fullfile(bsDB, bsProtocol, 'data', subject));
bsDataDir = fullfile(bsDB, bsProtocol, 'data', subject, d(end).name);
bsAnatDir = fullfile(bsDB, bsProtocol, 'anat', subject);

% 
if highres
    subFolder = 'highres';
    headmodelFileName = 'headmodel_surf_os_meg_02.mat';
else
    subFolder = 'lowres';
    headmodelFileName = 'headmodel_surf_os_meg.mat';
end

% create folder to save figures
saveDir  = fullfile(megRetDir, 'Quality_check', subject, 'BensonAtlasPred', subFolder);
if ~exist(saveDir, 'dir'); mkdir(saveDir); end


%% 1.Create combined hemi template that matches vertices in BS Gain matrix 

% Note: time consuming and computationally heavy if you want a high resolution template!
if (highres) && (~exist(fullfile(bsAnatDir, 'highres', 'benson14areas_overlay.mat')))
    interp_retinotopy(bsDB, fsdir, subject, subject, bsProtocol, highres)
end


if (~highres) && (~exist(fullfile(bsAnatDir, 'lowres', 'benson14areas_overlay.mat')))
    interp_retinotopy(bsDB, fsdir, subject, subject, bsProtocol)
end


%% 2. Load retinotopy templates

% Load V1-3 template with unitized phases in downsampled brainstorm format (computed by interp_retinotopy.m)
areas    = load(fullfile(bsAnatDir, subFolder,'benson14areas_overlay.mat')); % [1xNsources] Every value between [-3,3] is inside V1-3, zeros refer to outside of visual cortex. CHECK: Positive values represent LH (?) Negative values RH (?)
eccen    = load(fullfile(bsAnatDir, subFolder,'benson14eccen_overlay.mat')); % [1xNsources] Nonzero value represents vertex preferred eccentricity in degrees, zeros refer to outside of visual cortex
polarang = load(fullfile(bsAnatDir, subFolder,'benson14angle_overlay.mat')); % [1xNsources] Values represents vertex preferred polar angle in degrees (0 is the upper vertical meridian), zeros refer to outside of visual cortex
sigma    = load(fullfile(bsAnatDir, subFolder,'benson14sigma_overlay.mat')); % [1xNsources] Values represents vertex preferred polar angle in degrees (0 is the upper vertical meridian), zeros refer to outside of visual cortex

% get areas
area.v1 = (areas.sub_bs_areas==1);
area.v2 = (areas.sub_bs_areas==2);
area.v3 = (areas.sub_bs_areas==3);


% plot prf size by eccentricity
fn = fieldnames(area);
figure; hold all;

colorLines = hsv(length(fn));
for ii = 1:length(fn)
    theseVertices = area.(fn{ii});
    tmp = [eccen.sub_bs_eccen(theseVertices), sigma.sub_bs_sigma(theseVertices)];
    [val, ia, idx]  = unique(tmp(:,1), 'stable');
    binsize = sum(tmp(:,1)==unique(tmp(:,1))');
    
    mnEccenVsSize.(fn{ii}) = sort([val(:), accumarray( idx, tmp(:,2), [], @mean ) ], 1);
    seEccenVsSize.(fn{ii}) = sort([val(:), accumarray( idx, tmp(:,2), [], @std ) ./binsize ], 1);
    
    lo = mnEccenVsSize.(fn{ii})(:,2) - seEccenVsSize.(fn{ii})(:,2);
    hi = mnEccenVsSize.(fn{ii})(:,2) + seEccenVsSize.(fn{ii})(:,2);
    
    color = [0.5 0.5 0.5];
    err = patch([mnEccenVsSize.(fn{ii})(:,1); fliplr(mnEccenVsSize.(fn{ii})(:,1))], [lo; fliplr(hi)], color, 'FaceAlpha', 0.5);
%     plot(mnEccenVsSize.(fn{ii})(:,1), mnEccenVsSize.(fn{ii})(:,2), 'Color', colorLines(ii,:), 'LineWidth', 2);

end

%% 3. Get x,y coords from polar angle and eccen

areas.all = areas.sub_bs_areas>0;

% note, RH hemi should have negative polar angle values
theta = pi/180 * (90 - polarang.sub_bs_angle);
x = eccen.sub_bs_eccen .* cos(polarang.sub_bs_angle);
y = eccen.sub_bs_eccen .* sin(polarang.sub_bs_angle);

y0 = y(areas.all);
x0 = x(areas.all);

%% 4. Get RF with x,y,sigma and MEG stim
sigma = sigma.sub_bs_sigma(areas.all);

stim  = load(fullfile(megRetDir, 'Subject_sessions', subject, 'stimuli', 'meg', 'imported_stimulus', 'meg_stimulus.mat'));
stim = stim.meg_stim;

RF = rfGaussian2d(stim.X,stim.Y,sigma,sigma,false, x0, y0);

%% 5. Compute predicted BOLD response for every vertex in V1-V3, where pRF has unit volume
% Get volume height betas
beta = 1./sqrt(2*pi*sigma.^2);

% Predicted BOLD responses:
predBOLDV123 = bsxfun(@times, stim.im' * RF,beta');

predBOLD = zeros(size(areas.sub_bs_areas,1), size(predBOLDV123,1));
predBOLD(areas.all,:) = predBOLDV123';

clear RF x y
%% 6. Load Gain matrix and set dipole perpendicular to surface

% Only keep data sensors in Gain matrix
keep_sensors    = logical([ones(157,1); zeros(192-157,1)]); % Note: Figure out a more generic way to define keep_sensors

% Load headmodel from Brainstorm
headmodel = load(fullfile(bsDataDir, headmodelFileName));

% Get Gain matrix and truncate to not-nan sensors
G = headmodel.Gain(keep_sensors,:); % [Nsensors x 3*Nvertices]

% Contrained gain matrix
G_constrained = bst_gain_orient(G, headmodel.GridOrient); % [Nsensors x Nsources], equivalent to size BS pial cortex [1x15002]
clear G


%% 7. Load MEG data, get prediction, compare

% Get prediction
MEG.pred = G_constrained * predBOLD;

% Get data and conditions
tmp = load(fullfile(megRetDir, 'Subject_sessions', subject, 'data', 'meg', 'preproc', 'pp', 'epoched_data_hp_preproc_denoised.mat'));
MEG.data = tmp.data.data;
conditions = load(fullfile(megRetDir, 'Subject_sessions', subject, 'data', 'meg', 'preproc', 'pp', 'megStimConditions.mat'));
blinkIdx = conditions.triggers.stimConditions<20;

% Remove blink periods
MEG.predNoBlinks = MEG.pred;
MEG.predNoBlinks(:, ~blinkIdx(1:140)) = NaN;

% Get amplitudes at 10 Hz
MEG.F = fft(MEG.data);
flickerFreq = 10; % Hz
fs = 1000; % Hz
freqIdx = mprfFreq2Index(size(MEG.F,1), flickerFreq, fs);


chan = 20;
MEG.amps10 = nanmean(abs(MEG.F(freqIdx, :, :, chan)),3);

% Get modelfit data
d = dir(fullfile(megRetDir, 'Subject_sessions', subject, 'modeling', 'results', 'original_model', 'Run*', 'Results.mat'));
MEG.modelData = load(fullfile(d.folder, d.name));
MEG.datafit = MEG.modelData.results.orig_model.fit_data;

% Plot data and prediction
t = linspace(0, size(MEG.predNoBlinks,2).*1.1, size(MEG.predNoBlinks,2));

%%

fH = figure; set(gcf, 'Color', 'w', 'Position', [788, 798, 1146, 540]);

for sensor = 1:size(MEG.datafit,2)
    cla;
    
    plot(t, MEG.predNoBlinks(sensor,:), 'r', 'LineWidth', 2); hold on;
    % plot(t, MEG.amps10.*10^7, 'ko-')
    ylim([-1, 1].*10^-3); xlim([0,max(t)]);
    title(sprintf('Sensor %d - Predicted MEG response from Benson V1-V3 template', sensor))
    xlabel('time (s)'); ylabel('MEG response (AU)');
    set(gca, 'FontSize', 15, 'TickDir', 'out'); box off; 
    plot(t, MEG.datafit(:,sensor).*10^10, 'ko-')
    
    print(gcf, '-painters', '-depsc', fullfile(saveDir,sprintf('%d_bensonV123Pred', sensor)))
end
%% 
