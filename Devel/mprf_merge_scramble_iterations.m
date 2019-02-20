% mprf_merge_scramble_iterations

results_dir = fullfile(mprfRootPath,'data', 'subjectSession', 'wlsubj040', 'results', 'scramble_prfs');

d = dir(fullfile(results_dir, 'Run*'));

scrambled_corr = [];
orig_corr = [];


for ii = 1:length(d)
    r = load(fullfile(d(ii).folder, d(ii).name, 'Results.mat'));  
    
    scrambled_corr = cat(2, scrambled_corr, r.results.scrambled_corr);
    orig_corr = cat(2, orig_corr, r.results.orig_corr);

end

results = [];
results.scrambled_corr = scrambled_corr;
results.orig_corr = orig_corr;

model = r.model;

save(fullfile(results_dir, ['merged_' d(1).name]), 'results', 'model');