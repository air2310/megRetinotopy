sep_dir = '/Volumes/server/Projects/MEG/Retinotopy/Stimuli/Flicker_separate_ff';
save_dir = '/Volumes/server/Projects/MEG/Retinotopy/Stimuli/Flicker_frequency';

in_dir = dir(fullfile(sep_dir,'*.mat'));
fnames = {in_dir.name};

carriers = {'patterns','checkerboard','dartboard'};

flicker_freqs   = [30 20 15 12 10 6 5 4 ]; % Hz, i.e. images per second
n_perms = 30;
for n = 1:length(carriers);
    cur_car = carriers{n};
    
    cur_car_files = fnames(cellfind(strfind(fnames,cur_car)));
    
    if length(cur_car_files) == length(flicker_freqs)
    else
        error('Not all files found, or found too many')
    end
    
    tmp_stimuli = cell(size(cur_car_files));
    for nn = 1:length(cur_car_files)
        tmp_stimuli{nn} = load(fullfile(sep_dir,cur_car_files{nn}));
                
    end
    
    for nn = 1:n_perms
        cur_perm = randperm(length(cur_car_files));
        cur_stimuli = tmp_stimuli(cur_perm);
    
        fixSeq = [];
        trigSeq = [];
        diodeSeq = [];
        seq = [];
        seqtiming = [];
        
        for nnn = 1:length(cur_stimuli);
           
            fixSeq = [fixSeq; cur_stimuli{nnn}.stimulus.fixSeq];
            trigSeq = [trigSeq; cur_stimuli{nnn}.stimulus.trigSeq];
            diodeSeq = [diodeSeq; cur_stimuli{nnn}.stimulus.diodeSeq];
            seq  = [seq; cur_stimuli{nnn}.stimulus.seq];
            
            if nnn == 1
                to_add = 0;
            else
                to_add = max(seqtiming) + cur_stimuli{nnn}.stimulus.seqtiming(2);
            end
            
            seqtiming = [seqtiming; (cur_stimuli{nnn}.stimulus.seqtiming + to_add)];
            
            
        end
        
    stimulus = cur_stimuli{nnn}.stimulus;
    stimulus.fixSeq = fixSeq;
    stimulus.seqtiming = seqtiming;
    stimulus.trigSeq = trigSeq;
    stimulus.diodeSeq = diodeSeq;
    stimulus.seq = seq;
    
    save_path = fullfile(save_dir,['MEG_FF_' cur_car '_rand_perm_' num2str(nn)]);
    
    save(save_path,'stimulus');
        
    
    end
    
    
    
    
    
    
    
    
end
























