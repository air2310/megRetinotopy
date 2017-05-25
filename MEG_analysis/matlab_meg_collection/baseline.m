function data = baseline(datastack, pre, smpr)
%for each epoch in the datastack, for each channel in the data, takes the
%average of the values at that channel in that epoch for the period for
%pre data slices from the beginning of the epoch and subtracts that value
%from all the values for that channel for that epoch
data = zeros(length(datastack(:,1,1)),length(datastack(1,:,1)),length(datastack(1,1,:)));
pre = pre * smpr;
for i = 1:length(datastack(1,1,:))
    for j = 1:length(datastack(1,:,1))
        bl = mean(datastack(1:pre,j,i));
        data(:,j,i) = datastack(:,j,i) - bl;
    end;
end;
