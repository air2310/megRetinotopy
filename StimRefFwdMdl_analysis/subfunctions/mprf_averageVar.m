function averageDataToPlot = mprf_averageVar(ii,dataToPlot)
% function to compute the mean of the top 10 sensors in the back of the
% head

averageDataToPlot          = nanmean(dataToPlot(ii,:),1);

end
