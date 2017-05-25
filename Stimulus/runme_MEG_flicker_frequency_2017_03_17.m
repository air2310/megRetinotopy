function runme_MEG_flicker_frequency_2017_03_17(run, carrier)





tr = 1.5;
stimfile = sprintf('MEG_FF_%s_rand_perm_%d.mat',carrier,run);
fprintf('Loading: %s \n',stimfile);

cal = 'meg_lcd';
d   = loadDisplayParams(cal);
%         if any(size(d.gamma) == 1)
%             d.gamma = reshape(d.gamma,length(d.gamma)/3,3);
%         end


params = retCreateDefaultGUIParams;         % Some default parameters
% Check what these lines below do to the code...
params.modality         = 'meg';
params.prescanDuration  = 0;
params.interleaves      = NaN;
params.tr               = tr;
params.calibration      = cal;
params.framePeriod      = tr;
params.startScan        = 0;
params.motionSteps      = 2;
params.tempFreq         = 6/tr;
params.repetitions      = 1;
params.experiment       = 'Experiment From File';
params.period           = 12*params.tr;
params.numCycles        = 6;
params.skipSyncTests    = 1;
params.triggerKey       = '`';
params.fixation         = 'dot with grid';

params.loadMatrix = stimfile;
ret(params)

f = dir('~/Desktop/2017*.mat');
load(fullfile('~', 'Desktop', f(end).name));

% desired inter-stimulus duration
desired_timing = stimulus.seqtiming';
true_timing = response.flip - response.flip(1);
lag = true_timing - desired_timing;

figure;
plot(desired_timing,'b');
hold on;
plot(true_timing,'r');
legend('Desired timing','True timing');
ylabel('Seconds')
xlabel('Frame')

figure;
plot(lag)
ylabel('Seconds')
xlabel('Frame')
title('Lag (true timing - desiredtiming)');

figure;
plot(diff([0 lag]))
title('Lag increments (diff([0 lag]))')
ylabel('Seconds')
xlabel('Frame')


figure(101); clf

% desired inter-stimulus duration
plot(diff(stimulus.seqtiming));

% measured inter-stimulus duration
hold on; plot(diff(response.flip), 'r-');

ylim(median(diff(response.flip)) + [-.001 .001])
% frames between stimuli
frames = round(diff(response.flip) / (1/60));

% how many interstimulus frames differed from the median?
disp(sum(frames ~= median(frames)))


end

















