function eyeData = mprf_eyeAnalysis(subjID, dirPth, opt, modality, convertEDFFile)
% Function to load in eye tracking data from MEG or MRI session.
% Extract eye fixation and microsaccades for all blank and bar stimulus TRs
% or epochs.

%% 0.0 Define eyetracking analysis params
vThres      = 6; % Velocity threshold
msMinDur    = 6; % Microsaccade minimum duration (ms)
endTime     = 0; % number of messages to omit from end (if any)
blinkWinSec = [0.2 0.35]; % in seconds, window to define a blink (thus exclude from analysis)


%% 0.1 Convert EDF to mat file if requested
if convertEDFFile
    if strcmp(modality, 'MEG')
        filePth = dirPth.meg.eyePth;
    elseif strcmp(modality, 'MRI')
        filePth = dirPth.fmri.eyePth;
    end
    edffile = dir(fullfile(filePth, '*.edf'));
    for ii = 1:length(edffile)
        eyd = mglEyelinkEDFRead(fullfile(edffile(ii).folder, edffile(ii).name));
        
        save(sprintf(fullfile(filePth, '%s_eyelink.mat'),modality), sprintf('eyd%d',ii));
    end
end


%% 1. MEG analysis
if strcmp(modality, 'MEG')
    
    % Load eye link mat file and stimulus triggers
    load(fullfile(dirPth.meg.eyePth,sprintf('%s_eyelink.mat',modality)),'eyd');
    load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'megStimConditions.mat'), 'triggers');
    
    % Get start time in eyelink data
    for ii = 1:length(eyd.messages)
        if strfind(eyd.messages(ii).message, 'MEG Trigger')
            startTime = ii; break;
        end
    end
    
    % Define params, get xy data and messages
    timeLims    = [eyd.messages(startTime).time(1), eyd.messages(end-endTime).time];
    s           = msGetEyeData(eyd,timeLims,blinkWinSec);
    
    % Find start of run (20 is blink block)
    blinkTrigNr   = 20;
    blinkIndex    = find(triggers.stimConditions==blinkTrigNr);
    
    runStartIndex    = blinkIndex(1:2:end);
    diffNrOfEpochs   = diff(runStartIndex);
    nrOfEpochsPerRun = diffNrOfEpochs(1); % note, after the first blank epoch
    
    % Define other run details
    startEpoch          = opt.meg.epochStartEnd(1)*s.eyeInfo.smpRate;
    endEpoch            = opt.meg.epochStartEnd(2)*s.eyeInfo.smpRate;
    epochDur            = endEpoch-startEpoch; % should be 1100 ms with a sampling rate of 1000 Hz
    
    % Find the messages with run start
    blinkMessageIndex      = find(strcmp({eyd.messages.message}, sprintf('MEG Trigger: %d',blinkTrigNr)) == 1 ); % every run starts with message - MEG trigger: 20
    runStartMessageIndex   = blinkMessageIndex(1:2:end);
    
    % Get time points
    allTimePoints = [];
    xPos = []; yPos = [];
    xVel = []; yVel = [];
    for runIdx = 1:length(runStartMessageIndex)
        
        if mod(runIdx,6) % after 5 bar passes, 6th run there is eye movement (there was a break ?) - so that was excluded
            
            startOfRun = runStartMessageIndex(runIdx);
            
            startOfRun = startOfRun+2; % skip the two blink blocks
            
            epochsToInclude  = startOfRun:(startOfRun+nrOfEpochsPerRun-3); % 25 epochs (27 excluding: 3 blanks and 22 bar passes minus 2 blinks)
            
            epochsToInclude_starttime = [eyd.messages(epochsToInclude).time] - s.timeRaw(1); % time is calculated with respect to a reference time
            
            thisRunTimePoints = epochsToInclude_starttime'+(startEpoch:endEpoch-1); % there are 1100 points between consecutive epochsToInclude_time values, but each epoch was 1300 ms long
            
            allTimePoints = cat(1, allTimePoints, thisRunTimePoints); % total should be (95x25)x1100, or (5 sweeps x 19runs) by 25 epochs per sweep by 1100 timepoints
            
            % Get eye position and velocity for selected timepoints in epochs
            temp_xPos = s.xyPos(thisRunTimePoints(:),1);
            xPos = cat(1, xPos, reshape(temp_xPos, size(thisRunTimePoints,1), size(thisRunTimePoints,2)));
            
            temp_yPos = s.xyPos(thisRunTimePoints(:),2);
            yPos = cat(1, yPos, reshape(temp_yPos, size(thisRunTimePoints,1), size(thisRunTimePoints,2)));
            
            temp_xVel = s.xyVel(thisRunTimePoints(:),1);
            xVel = cat(1, xVel, reshape(temp_xVel, size(thisRunTimePoints,1), size(thisRunTimePoints,2)));
            
            temp_yVel = s.xyVel(thisRunTimePoints(:),2);
            yVel = cat(1, yVel, reshape(temp_yVel, size(thisRunTimePoints,1), size(thisRunTimePoints,2)));
            
        end
    end
    
    %% 2. MRI analysis
elseif strcmp(modality, 'MRI')
    
    % Load eye link mat file and stimulus triggers
    eydData = load(fullfile(dirPth.fmri.eyePth,sprintf('%s_eyelink.mat',modality)));
    load(fullfile(dirPth.fmri.stimFilePth, 'MRI_retinotopy_stimulus_run_1.mat'), 'stimulus');
    
    % Define params, get xy data and messages
    endTime     = 0; % number of messages to omit from end (if any)
    blinkWinSec = [0.2 0.35];
    
    % Define other run details
    stimRad     = 10;   % deg
    TR          = 1500;  % msec
    nTRs        = 244;  % nr volumes
    epochTimes  = linspace(0, (nTRs-1)*TR, nTRs);
    runs        = 6;
    
    % Get start time in eyelink data
    s = cell(1,runs);
    xPos = []; yPos = [];
    xVel = []; yVel = [];
    
    for fn = 1:runs
        
        eyd       = eydData.(sprintf('eyd%d',fn));
        
        if ~isempty(eyd)
            
            timeLims  = [eyd.messages(end).time, eyd.gaze.time(end)];
            s{fn}     = msGetEyeData(eyd,timeLims,blinkWinSec);
            
            epochSampleRate = (TR/1000) * s{fn}.eyeInfo.smpRate; % sec
            
            % Pad data with zero's if necessary
            paddingLength = (nTRs*epochSampleRate) - length(s{fn}.xyPos(:,1));
            
            % Get eye position and velocity for selected timepoints in TR
            if paddingLength > 0
                temp_xPos = [s{fn}.xyPos(:,1); zeros(paddingLength,1)];
                xPos = cat(1, xPos, reshape(temp_xPos, epochSampleRate, nTRs));
                
                temp_yPos = [s{fn}.xyPos(:,2); zeros(paddingLength,1)];
                yPos = cat(1, yPos, reshape(temp_yPos, epochSampleRate, nTRs));
                
                temp_xVel = [s{fn}.xyVel(:,1); zeros(paddingLength,1)];
                xVel = cat(1, xVel, reshape(temp_xVel, epochSampleRate, nTRs));
                
                temp_yVel = [s{fn}.xyVel(:,2); zeros(paddingLength,1)];
                yVel = cat(1, yVel, reshape(temp_yVel, epochSampleRate, nTRs));
                
            else paddingLength < 0
                temp_xPos = s{fn}.xyPos(1:end+paddingLength,1);
                xPos = cat(1, xPos, reshape(temp_xPos, epochSampleRate, nTRs));
                
                temp_yPos = s{fn}.xyPos(1:end+paddingLength,2);
                yPos = cat(1, yPos, reshape(temp_yPos, epochSampleRate, nTRs));
                
                temp_xVel = s{fn}.xyVel(1:end+paddingLength,1);
                xVel = cat(1, xVel, reshape(temp_xVel, epochSampleRate, nTRs));
                
                temp_yVel = s{fn}.xyVel(1:end+paddingLength,2);
                yVel = cat(1, yVel, reshape(temp_yVel, epochSampleRate, nTRs));
            end
        end
    end
    
end


% Concatenate all conditions for XY position and XY velocity
xyPos = cat(3,xPos, yPos);
xyVel = cat(3,xVel, yVel);

if strcmp(modality, 'MRI') % put TRs as first dim (equivalent to epochs in MEG)
    xyPos = permute(xyPos, [2,1,3]);
    xyVel = permute(xyVel, [2,1,3]);
end


% Get microsaccades
epochsWithMS = [];
for epoch = 1:size(xyVel,1)
    
    theseXYPos = squeeze(xyPos(epoch,:,:));
    theseXYVel = squeeze(xyVel(epoch,:,:));
    
    % Define microsaccades
    [sacRaw,radius] = microsacc(theseXYPos,theseXYVel,vThres,msMinDur);
    
    % Remove the ones that occurr closely together (overshoot)
    numSacs = size(sacRaw,1);
    minInterSamples = ceil(0.01*s{1}.eyeInfo.smpRate);
    
    if numSacs ~= 0
        interSac = sacRaw(2:end,1)- sacRaw(1:end-1,2);
        sac = sacRaw([1; find(interSac > minInterSamples)+1],:);
        fprintf('%d rejected for close spacing\n', numSacs - size(sac,1));
        fprintf('%d saccades detected\n', size(sac,1));
        
        % Saved detected saccades into variable called 'sAll'
        sAll.sacsRaw(epoch)          = {sacRaw};
        sAll.sacs(epoch)             = {sac};
        sAll.sacDetectRadius(epoch)  = {radius};
        sAll.eyeInfo.vThres(epoch,:)   = vThres;
        sAll.eyeInfo.msMinDur(epoch,:) = msMinDur;
        sAll.numSacs(epoch)            = numSacs;
        
        
        fprintf('number of saccades detected: %d, detectRadius: [%0.3f %0.3f]\n', ...
            size(sAll.sacs{epoch},1), sAll.sacDetectRadius{epoch}(1), sAll.sacDetectRadius{epoch}(2));
        
        % Make a binary microsaccade vector where each 1 is the start
        % of a ms in an epoch (1 s)
        for msNr = 1:numSacs
            msVec(epoch,sacRaw(msNr,1)) = 1;
        end
        
        epochsWithMS = [epochsWithMS epoch];
        
        % Get circular distribution of ms
        [rho{epoch},theta{epoch}] = msGetThetaRho(squeeze(xyPos(epoch,:,:)),sAll.sacs{epoch});
        
    end
end

eyeData = struct();
eyeData.xPos    = xPos;     % deg per time point
eyeData.yPos    = yPos;     % deg per time point
eyeData.xVel    = xVel;     % deg per time point
eyeData.yVel    = yVel;     % deg per time point
eyeData.xyPos   = xyPos;    % deg per time point
eyeData.xyVel   = xyVel;    % deg per time point
eyeData.s       = s;        % eye data struct
eyeData.sAll    = sAll;     % struct with all ms data
eyeData.msVec   = msVec;    % vector, ms counts per epoch
eyeData.rho     = rho;      % ms count, for circular distribution
eyeData.theta   = theta;    % ms orientation, for circular distribution
eyeData.epochSampleRate = epochSampleRate; % sec

