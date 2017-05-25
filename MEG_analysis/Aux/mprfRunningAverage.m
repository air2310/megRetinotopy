function av = mprfRunningAverage(cur_val, new_val, n)


av = (cur_val .* ((n-1)./n)) + (new_val .* (1-(n-1)/n));




end