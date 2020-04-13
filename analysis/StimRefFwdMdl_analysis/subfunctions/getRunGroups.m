function runGroup = getRunGroups(nRuns)

% %If randomized
% tmp = randperm(nRuns);
% runGroup{1} = tmp(1:9);
% runGroup{2} = tmp(10:nRuns);

runGroup{1} = 1:2:nRuns;
runGroup{2} = 2:2:nRuns;

end