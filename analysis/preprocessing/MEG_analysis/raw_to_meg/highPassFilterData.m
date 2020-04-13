function filteredData = highPassFilterData(dataIn, fParams, verbose)

% Create design with high pass and band pass filter
h = fdesign.highpass(fParams.fStop, fParams.fPass, fParams.aStop, fParams.aPass, fParams.fs);
Hd = design(h,'butter','matchexactly','stopband');

% Visualize
% if verbose
%     fvtool(Hd);
% end

origDim = size(dataIn);

fprintf('Filtering data from %d channels:\n',origDim(3))
for n = 1:size(dataIn,3)
    fprintf('%d.',n);

    tmpData = dataIn(:,:,n);
    tmpData = tmpData(:);
    nanIdx = isnan(tmpData);
    
    dataOut = nan(size(tmpData));
    dataOut(~nanIdx) = filter(Hd,tmpData(~nanIdx));

    filteredData(:,:,n) = reshape(dataOut,[origDim(1:2) 1]);

end

fprintf('\nDone.\n');

end









