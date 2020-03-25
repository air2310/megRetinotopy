function splitHalfAmpReliability = mprf_ComputeSplitHalfReliability(dirPth, opt, megData,freqIdx)

% Check dimensions of MEG data
[~, ~, nRuns, nSensors] = size(megData);
nIter = 1000;

% Pre-allocate space
allVarExpl    =  NaN(2, nIter, nSensors);
allVarExpl_mn =  NaN(nIter, nSensors);

for n = 1:nIter
    
    tmp = randperm(nRuns);
    
    % Split runs in half
    r1 = tmp(1:9);
    r2 = tmp(10:nRuns);
    
    for s = 1:nSensors
        
        if opt.meg.useCoherentSpectrum
            
            % Average runs within half
            meanTimeseries1  = nanmean(megData(:,:,r1,s),3);
            meanTimeseries2  = nanmean(megData(:,:,r2,s),3);
            
            % Get amplitudes
            F1   = fft(meanTimeseries1);
            amp1_mn = abs(F1(freqIdx,:))/size(meanTimeseries1,1)*2;
            
            F2   = fft(meanTimeseries2);
            amp2_mn = abs(F2(freqIdx,:))/size(meanTimeseries2,1)*2;
            

        else
            % First to fft
            F   = fft(megData);
            amps = abs(F(freqIdx,:))/size(megData,1)*2;
            
            % then select runs
            amp1 = amps(:,r1,s)';
            amp2 = amps(:,r2,s)';
            
            % average them
            amp1_mn = nanmean(amp1);
            amp2_mn = nanmean(amp2);
            
             % Remove nans
            amp1_mn = amp1_mn(~isnan(amp1_mn));
            amp2_mn = amp2_mn(~isnan(amp2_mn));
        end
        
        CoD = @(data,prediction) (1 - (var(data - prediction, 'omitnan') ./ var(data, 'omitnan')));
        
        varexpl1 = CoD(amp1_mn'.*10^14, amp2_mn'.*10^14);
        varexpl2 = CoD(amp2_mn'.*10^14, amp1_mn'.*10^14);
%  
%         allVarExpl(1,n,s) = varexpl1;
%         allVarExpl(2,n,s) = varexpl2;
        
        allVarExpl_mn(n,s) = nanmean([varexpl1,varexpl2]);
            
    end
end

splitHalfAmpReliability = nanmedian(allVarExpl_mn,1);

save(fullfile(dirPth.model.saveFigPth, opt.subfolder, 'refphase', 'splitHalfAmpReliability1000.mat'), 'splitHalfAmpReliability');

return