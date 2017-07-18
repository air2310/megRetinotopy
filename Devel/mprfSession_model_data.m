



pred = mprf__load_model_predictions;
[fname, fpath] = uigetfile('*.mat', 'Select raw data to model');
% 
tmp = load(fullfile(fpath, fname));
var_name = fieldnames(tmp);
epoched_data = tmp.(var_name{1});
clear tmp

periods.blank = [3:5 30:32 57:59 84:86 111:113 138:140];
periods.blink = [1 2 28 29 55 56 82 83 109 110 136 137];
periods.stim = setdiff(1:140,[periods.blink periods.blank]);

%mprf__coranal_on_meg_data

sz = size(epoched_data.data);
n_time = sz(1);
n_bars = sz(2);
n_reps = sz(3);
n_chan = sz(4);

fprintf('Processing %d stimulus periods:\n',sz(2));

opts.measure = 'stim_amp';
opts.do_av = true;
opts.do_std = true;
opts.do_ste = true;

samp_rate = 1000;
stim_freq = 10;
opts.stim_idx = round(mprfFreq2Index(n_time,stim_freq,samp_rate));

for n = 1:n_bars
    fprintf('%d.',n)
    
    for nn = 1:n_chan
        cur_data = squeeze(epoched_data.data(:,n,:,nn));
        
        tmp = mprf__compute_meg_measure(cur_data,opts);
        
        measure.raw(:,n,nn) = tmp.raw;
        
        if opts.do_av
            if n == 1 && nn ==1
                measure.av = nan(n_bars, n_chan);
            end
            measure.av(n,nn) = tmp.av;
        end
        
        
        if opts.do_std
            
            if n == 1 && nn ==1
                measure.std = nan(n_bars, n_chan);
            end
            measure.std(n,nn) = tmp.std;
            
            
        end
        
        
        if opts.do_ste
            
            if n == 1 && nn ==1
                measure.ste = nan(n_bars, n_chan);
            end
            measure.ste(n,nn) = tmp.ste;
            
            
        end
        
        
    end
end
fprintf('\n')






