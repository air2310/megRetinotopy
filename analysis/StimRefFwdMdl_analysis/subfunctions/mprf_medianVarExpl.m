function medianDataToPlot = mprf_medianVarExpl(ii,dataToPlot)
% function to compute the median of the top 10 sensors

medianDataToPlot          = nanmedian(dataToPlot(ii,:),1);

end