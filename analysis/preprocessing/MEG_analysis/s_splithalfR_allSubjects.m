subjIDs = {'wlsubj004', 'wlsubj039', 'wlsubj040', 'wlsubj058', ...
    'wlsubj068', 'wlsubj070', 'wlsubj081', 'wlsubj106', ...
    'wlsubj109', 'wlsubj111'};

% Set options if not defined (see getOpts function for more options)
opt = getOpts('saveFig', true,'verbose', true, 'fullSizeMesh', true, ...
    'perturbOrigPRFs', false, 'addOffsetParam', false, ...
    'refitGainParam', false);

for s = 1:length(subjIDs)
    
    subjID =  subjIDs{s};
    
    % Load paths with data files for this subject
    dirPth = loadPaths(subjID);
    
    fprintf('(%s): Starting analysis of subject %s, using %s\n', mfilename, subjID, regexprep(opt.subfolder,'/',' '));
    
    % Load MEG data
    load(fullfile(dirPth.meg.processedDataPth, 'allEpochs', 'epoched_data_hp_preproc_denoised.mat'), 'data');
    megData = data.data;
    
    % Check dimensions of MEG data
    [nTimepoints, nEpochs, nRuns, nSensors] = size(megData);
    
    % Get freq index
    freqIdx = mprfFreq2Index(nTimepoints, opt.meg.flickerFreq, opt.meg.fs);
    
    [splitHalfAmpReliability(s,:),splitHalfAmpCorrelation(s,:)] = mprf_ComputeSplitHalfReliability(dirPth, opt, megData,freqIdx);
end

