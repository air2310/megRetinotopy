function trigger = meg_fix_triggers_wlsubj058(ts, triggerChan)

%% 0. Part where there is a normal trigger sequence

ts3 = ts(triggerChan(1:5),:);
ts3(:,1:3688049)=1;
trigger3 = meg_fix_triggers(ts3');

trig3_ind = find(trigger3);

condsOneRun = trigger3(trig3_ind(1:140));



%% 1. Part where every trigger has the same value

ts1 = ts(triggerChan(1),:);
ts1(1,1:110236) = 0;
ts1(1,3376547:end) = 0;

md = median(ts1(1,110236:3376546));
ts1(ts1 > 5*md) = NaN;

 
% rescale to [0 1]
ts1 = ts1 - min(ts1(:));
ts1 = ts1 / max(ts1(:));

% check whether triggers are indicated by a low value or a high value
if round(mean(ts1(:))) == 0, trigger_is_high = true; 
else                     trigger_is_high = false; end 

% if triggers are indicated by a low value, then invert the signal
if ~trigger_is_high, ts1 = 1 - ts1; end

% threshold to binarize from analog recording 
ts1 = ts1 > 0.1;
% differentiate to isolate trigger onsets, and pad to preserve length
ts1 = padarray(diff(ts1'), [1 0], 0, 'post');

% rectify 
ts1(ts1<0) = 0;

% mark the samples as 1 when any trigger channel signaled 
any_trigger      = sum(ts1,2) > 0;

% list the sample numbers (time ponts) when any trigger channel signalled
any_trigger_inds = find(any_trigger);

% check whether 2 triggers signalled in an impossibly short time period
% (indicating that they were probably simultanuous, but slighly misaligned)
triggers_that_are_too_close = diff(any_trigger_inds) < 10;

time_points_that_are_too_close = any_trigger_inds(triggers_that_are_too_close);

% Issue: most of the above indexes return all zeros, check on that 
% Solution: because of padding, the triggers are present -1 or -2 time
% points before the the time points listed in bad_time_points

% Realign bad triggers
%   We find any case where two trigger channel onsets were within 2
%   samples, and we re-align them to the first sample.
for ii=1:length(time_points_that_are_too_close)
    ts1(time_points_that_are_too_close(ii),:) = sum(ts1(time_points_that_are_too_close(ii)+[-2:2],:));
    ts1(time_points_that_are_too_close(ii) + [-2 -1 1 2],:) = 0;
end

% now check that there are no time points that are too close
any_trigger                 = sum(ts1,2) > 0;
any_trigger_inds            = find(any_trigger);
triggers_that_are_too_close = diff(any_trigger_inds) < 10;
if sum(triggers_that_are_too_close) > 0
    ts1(time_points_that_are_too_close) = 0;
end

any_trigger_inds2  = find(ts1);
toRemove = any_trigger_inds2((find(diff(find(ts1))<1250)+1));
toRemove = [toRemove', 1023877, 2327021, 331921, 548872, 764238, 978468, 1013709, 1228852, 1438635, 1650593, 1861832, 2080989, 2315871, 2527878, 2749867, 2964736, 3178934];

ts1(toRemove) = 0;

% convert binary matrix into base 10 vector
trigger1 = ts1 * 2.^(0:size(ts1,2)-1)';

conditions1= repmat(condsOneRun, [1 15]);
trigger1(find(trigger1)) = conditions1(:);

%% 2. Part where there is no trigger, so use Photodiode

pd_chan = 192;

% Get peaks of photodiode response
dataPD  = diff(ts(pd_chan,:));

dataPD(1:3450486)=0;
dataPD(3634693:end)=0;

dataPD = padarray(dataPD', [1 0], 0, 'post');
dataPD = dataPD';

threshPD     = (max(dataPD) - min(dataPD)) / 10;
[pd.value_white, pd.ind_white] = findpeaks([0, -dataPD],'MINPEAKHEIGHT',threshPD);
[pd.value_black, pd.ind_black] = findpeaks([0, dataPD],'MINPEAKHEIGHT',threshPD);

pd.ind = sort([pd.ind_white, pd.ind_black]);

[val, idx] = find(diff(pd.ind_black) >1290);

trig_ind = pd.ind_black(idx);

trig_ind = [trig_ind, pd.ind_black(idx(end))+1300];

trigger2 = zeros(size(dataPD));
trigger2(trig_ind) = condsOneRun;

%% Combine them all
trigger = [trigger1 + trigger2' + trigger3];

end


