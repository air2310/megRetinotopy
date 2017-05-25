function [flip_time, stim_time, init, time0] = mprfGetStimAndFlipTime(stimulus, response, params, time0)

if isempty(time0)
    time0 = response.flip(1);
end

flip_time.first_flip = response.flip(1);
flip_time.last_flip = response.flip(end);

stim_time.seq_times = ((stimulus.seqtiming + flip_time.first_flip - time0)*1000);% The time for which a flip was requested
% The code presents all frames in the sequence and does not
% drop any of them, or does not speed the stimulus presentation
% up to keep up with the timing. Therefore, missing the VBL
% deadline will cause the flips to be delayed relative to the
% stimulus sequence timing.

flip_time.flip_times = ((response.flip - time0)...
    .* 1000);% The times at which the screen flipped

flip_time.trigger_times = round((response.flip(stimulus.trigSeq > 0)...
    - time0).*1000); % Times at which a flip corresponding to a trigger happened

flip_time.trigger_times_02 = zeros(size(stim_time.seq_times));

if isfield(response,'trig')
    flip_time.trigger_times_02(1:length(response.trig))...
        = response.trig'; % The trigger id at its corresponding flip
else
    warning('No .trig field in response variable');
end

stim_time.trigger_times = round((stimulus.seqtiming(stimulus.trigSeq > 0)...
    +flip_time.first_flip - time0)*1000); % Times at which a trigger was requested

stim_time.trigger_idx = stimulus.trigSeq(stimulus.trigSeq > 0);

init.flip_times = round((stimulus.flashTimes.flip - time0) * 1000);
init.seq = params.display.initstim;

end









