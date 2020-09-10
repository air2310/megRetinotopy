%% makeFigureXX_eyeMovements_xyPosition
% Script to plot horizontal and vertical eye movement for MEG scan run for
% MEG retinotopy project.
%
% Note: This function relies on code from noisepoolPCADenoise repository

% For MRI data choose from wlsubj039, wlsubj058, wlsubj068, wlsubj070, wlsubj109, wlsubj111
% For MEG data choose from wlsubj039, wlsubj058, wlsubj068, wlsubj070, wlsubj081, wlsubj106, wlsubj109, wlsubj111
modality = 'MEG';       % choose 'MEG' or 'MRI' (Note: subjects 004, 040 have neither eye datasets, 081 has no MRI eye data)
subjID   = 'wlsubj039'; 
equateSamples = false;  % equate number of samples for MEG to match MRI 

% Load paths with data files for this subject
dirPth = loadPaths(subjID);

% Go back to root
cd(mprf_rootPath)

% Set options
opt = getOpts('saveFig',1,'verbose',1); % see getOpts function for more options

% convert raw EDF file or use saved out mat file?
convertEDFFile =  false; 

eyeData = mprf_eyeAnalysis(subjID, dirPth, opt, modality, convertEDFFile);

% If equating the number of samples for MRI and MEG, then randomly select
% 1000 epochs:
%   MRI => 244 epochs x 4500 samples (500 Hz sample rate) = 1,098,000 samples
%   MEG => 2375 epochs x 1100 samples (1000 Hz sample rate) = 2,612,500
%   samples
% Thus ~1,100,000 MRI samples / 1100 samples in one MEG epoch = ~1000 epochs
if strcmp(modality, 'MEG') && equateSamples
    rn_vector = randi(size(eyeData.msVec,1),1,1000);
    eyeData.xyPos = eyeData.xyPos(rn_vector,:,:);
    eyeData.xyVel = eyeData.xyVel(rn_vector,:,:);
    eyeData.xPos  = eyeData.xPos(rn_vector,:);
    eyeData.yPos  = eyeData.yPos(rn_vector,:);
    eyeData.msVec  = eyeData.msVec(rn_vector,:);
    eyeData.theta  = eyeData.theta(rn_vector);
    eyeData.rho  = eyeData.rho(rn_vector);

end

% Get mean, median and covariance matrix
traceMean   = nanmean([eyeData.xPos(:), eyeData.yPos(:)]);
traceMedian = nanmedian([eyeData.xPos(:), eyeData.yPos(:)]);
traceCov    = cov([eyeData.xPos(~isnan(eyeData.xPos(:))), eyeData.yPos(~isnan(eyeData.yPos(:)))]);

fprintf('Diagonal of covariance matrix %1.2f\t%1.2f\n',diag(traceCov))

%% -------------------------------------
% ------------- Fixations --------------
% --------------------------------------

% plot horizontal vs vertical eyemovement positions
fH1 = figure(1); clf;
set(gcf,'Color', 'w', 'Position',[31, 252, 1067, 996]);  hold on;
% Plot eye traces
plot(eyeData.xPos,eyeData.yPos, '.','MarkerSize',1, 'Color', [.7 .7 .7]);
% Plot 95% confidence ellipse
ax = error_ellipse(traceCov,traceMean,'conf',0.95,'color','k');
% Plot median
plot(traceMean(1),traceMean(2), 'r+', 'markersize',5);

% Scale axes
stimRad = 10; %deg
xlim([-stimRad stimRad]);
ylim([-stimRad stimRad]);
xlabel('X Position (deg of visual angle)');
ylabel('Y Position (deg of visual angle)');
grid on; box off;
axis square;
title(sprintf('Subject %s - Fixation %s', subjID, modality))
set(gca,'FontUnits','centimeters', 'FontSize', 0.5,'TickDir','out','LineWidth',2);


% save the figure
if opt.saveFig
    [pth, ~] = fileparts(dirPth.model.saveFigPth);
    saveDir = fullfile(pth, 'finalfig', 'figureXX_EyeFixation');
    if ~exist(saveDir, 'dir')
        mkdir(saveDir);
    end
    fprintf('\n(%s): Saving eye movements figure  in %s\n',mfilename, saveDir);
    
    if strcmp(modality, 'MEG') && equateSamples
        fname = sprintf('figXX_eyeMovements_%s_%s_equateSamples', subjID,lower(modality));
    else
        fname = sprintf('figXX_eyeMovements_%s_%s', subjID,lower(modality));
    end

    print(fH1, fullfile(saveDir, fname), '-dpng');
    print(fH1, fullfile(saveDir, fname), '-depsc');
    
end

%% -------------------------------------
% ------------- MS rates ---------------
% --------------------------------------

fH2 = figure(100); clf; set(gcf, 'Color', 'w'); hold all;

nboot    = 1000;
maxbin   = median(sum(eyeData.msVec,2))+5;
nbins    = linspace(0,maxbin, 61);
epochnum = 1:size(eyeData.msVec,1);
bootstat = bootstrp(nboot, @(x) nanmean(sum(eyeData.msVec(x, :),2)), epochnum);

mdBootsStat = median(bootstat);
y = hist(bootstat,nbins);
plot(nbins, y, 'Color', 'k','LineWidth',2);
plot([mdBootsStat mdBootsStat],[0 500],'--','Color', 'k');

fprintf('Median of bootstrapped ms rate %1.2f\n',mdBootsStat)


xlabel('MS rate/s')
ylabel('Frequency')
makeprettyaxes(gca, 18,18);
title(sprintf('MS Distributions subject %s %s', subjID, modality));

if opt.saveFig
    if strcmp(modality, 'MEG') && equateSamples
        fname = sprintf('figXX_msRateHistogram_%s_%s_equateSamples', subjID,lower(modality));
    else
        fname = sprintf('figXX_msRateHistogram_%s_%s', subjID,lower(modality));
    end
    
    print(fH2, fullfile(saveDir, fname), '-dpng');
    print(fH2, fullfile(saveDir, fname), '-depsc');
end



%% -------------------------------------
% ------ Circular MS histograms --------
% --------------------------------------

idx = cellfun(@isempty, eyeData.theta);
theta_toplot = eyeData.theta(~idx);
rho_toplot   = eyeData.rho(~idx);

thetas = [];
for ii = 1:numel(theta_toplot)
    thetas = vertcat(thetas, theta_toplot{ii});
end

rhos = [];
for ii = 1:numel(rho_toplot)
    rhos = vertcat(rhos, rho_toplot{ii});
end

% plot individual microsaccade vectors
% h=polar2(thetas,rhos,radialLimit,'o');
% set(h,'markersize',2);
% set(h,'Color','k');
% title(sprintf('Microsaccade vectors S%s',subjID))

% plot circular distribution
nbins = length(0:10:360);
fH3 = figure; set(gcf, 'Color', 'w'); polarhistogram(thetas, nbins, 'FaceColor', 'k', 'FaceAlpha', 0.4)

set(gca, 'LineWidth', 2, 'FontSize',14)
title(sprintf('Angular distribution MS - subject %s %s',subjID, modality))


if opt.saveFig
    if strcmp(modality, 'MEG') && equateSamples
        fname = sprintf('figXX_msHistogramAngular_%s_%s_equateSamples', subjID,lower(modality));
    else
        fname = sprintf('figXX_msHistogramAngular_%s_%s', subjID,lower(modality));
    end
    print(fH3, fullfile(saveDir, fname), '-dpng');
    print(fH3, fullfile(saveDir, fname), '-depsc');
end



