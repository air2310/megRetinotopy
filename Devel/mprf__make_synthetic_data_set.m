function mprf__make_synthetic_data_set(syn, meg_resp, cur_time, channels)

main_dir = mprf__get_directory('main_dir');
meg_ddir = mprf__get_directory('meg_data');
param_dir = fullfile(main_dir,mprf__get_directory('meg_stim_params'));

template_file = fullfile(main_dir, meg_ddir, syn.file);
hdr = ft_read_header(template_file);

fprintf('Getting timing info...\n')
timing = mprfGetTriggers(template_file, param_dir, channels.triggers, channels.diode);

try
    triggers = [timing.trigger.trigger2flip(:,1) timing.trigger.idx(:)];
    
catch
    triggers = [timing.trigger.channel(:,1) timing.trigger.idx(:)];
    
end

fprintf('Done.\n');


syn_dir = fullfile(main_dir,mprf__get_directory('syn_data'));
if ~exist(syn_dir, 'dir')
    mkdir(syn_dir)
end

save_name = fullfile(syn_dir, ['syn_data_modelling_run_' cur_time '.sqd']);
save(fullfile(syn_dir, ['trig_info_modelling_run_' cur_time '.mat']),'triggers','timing');

n = 1300; % time points per epoch
t = (1:n)/1000;
stim_locked_signal = sin(t*2*pi*10);

noise_s = .2*10^4.75;
sl_s    = 5*10^4.75;
bb_s    = 1*10^4.75;


% for generation of pink noise
f       = -n/2:n/2-1;
alpha   = 1;
cpfov   = ifftshift(f);
amp     = 1./(abs(cpfov).^alpha);
amp(1)  = 0;


num_chan = 157;
num_runs = 19;
num_stims = size(meg_resp,1);

data157 = randn(num_chan, hdr.nSamples)*10^2.5;


for ii = 1:num_runs
    fprintf('.');
    for jj = 1:num_stims
        
        bb    = real(ifft(fft(randn(n,num_chan),[],1) .* repmat(amp',[1 num_chan]),[],1));
        noise = real(ifft(fft(randn(n,num_chan),[],1) .* repmat(amp',[1 num_chan]),[],1));
        
        idx = (0:1299)+ triggers((ii-1)*140+jj,1);
        sl      = meg_resp(jj,:)' * (stim_locked_signal .* sl_s);
        bb      = bsxfun(@times, bb, meg_resp(jj,:))' *  bb_s;
        noise   = noise' * noise_s;
        
        data157(:,idx) = sl + bb + noise;
    end
end

fprintf('\n')


fprintf('Copying data file...\n')
system(sprintf('cp %s %s',template_file,save_name))

fprintf('Done.\n')

max_chan_per_chunk = 60;
chunk_idx = 1:max_chan_per_chunk*ceil(num_chan / max_chan_per_chunk);
chunk_idx = reshape(chunk_idx,max_chan_per_chunk,[]);

for n = 1:size(chunk_idx,2);
    cur_idx = chunk_idx(:,n);
    
    if max(cur_idx) > size(data157,1)
        cur_idx = cur_idx(cur_idx <= size(data157,1));
        
        
    end
    
    fprintf('Writing channels:\n')
    fprintf('%d : %d\n',cur_idx(1), cur_idx(end))
    sqdwrite(save_name, save_name ,...
        'Channels', channels.data(cur_idx), ...
        'Data',data157(cur_idx,:)',...
        'Action','Overwrite',...
        'Samples',[1 size(data157,2)]);
end





end











