function mprfSyntheticData


project_dir = mprfRootPath;
project_dir = project_dir(1:strfind(project_dir, '/Code/mprfSession'));

cur_dir = pwd;
cd(project_dir);

[pred_fname, pred_dir] = uigetfile('*','Please select prediction file to use');
pred_file = fullfile(pred_dir,pred_fname);

main_sub_dir = 'Subject_sessions';

tmp_sd = strfind(pred_dir,main_sub_dir) + length(main_sub_dir)+1;
tmp_fb = find(pred_dir(tmp_sd:end) == '/',1,'first')-1;
sub_dir = pred_dir(1:tmp_sd + tmp_fb);

cd(sub_dir);

load(fullfile(pwd,'mprfSESSION.mat'));
global mprfSESSION


% Still hard coded for convenience. May what to ask user which data to use
% as .sqd template, replacing data_dir and raw_file. Also, may want to
% select param file directory manually. Additionally, make the
% mprfPreprocessRetinotopyData take in the timing file generated here...

%wl_subj040:
%data_dir = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wl_subj040/wl_subj040_20170406';   % Subject's data directory
%raw_file = 'R1151_Retinotopy_04.06.17.sqd';                             % Raw file
%source_dir = fullfile(data_dir,'Raw');
%param_dir = fullfile(data_dir,'Stimulus/Param_files');                  % Parameter files

%wl_subj_004:
data_dir = '/Volumes/server/Projects/MEG/Retinotopy/Data/MEG/wl_subj_004';   % Subject's data directory
raw_file = 'R0774_RETMEG_Comb_5.12.17.sqd';                             % Raw file
source_dir = fullfile(data_dir,'Raw','Export_02');
param_dir = fullfile(data_dir,'Stimulus/Param_files');                  % Parameter files



% Channel identity:
trig_chan = 160 : 167;
diode_chan = 191;
data_chan = 0:156;
raw_meg_file = fullfile(source_dir,raw_file);
hdr = ft_read_header(raw_meg_file);

fprintf('Getting timing info...\n')
timing = mprfGetTriggers(raw_meg_file, param_dir, trig_chan, diode_chan);

try
    triggers = [timing.trigger.trigger2flip(:,1) timing.trigger.idx(:)];
    
catch
    
    triggers = [timing.trigger.channel(:,1) timing.trigger.idx(:)];
    
end

fprintf('Done.\n');

pred = load(pred_file);

n = 1300; % time points per epoch
t = (1:n)/1000;
stim_locked_signal = sin(t*2*pi*10);

num_chan = 157;
num_runs = 19;
num_stims = size(pred.gain.M_02,1);

%data = sqdread(raw_meg_file,'Channels', 157:191);

% make simulation
data157 = randn(num_chan, hdr.nSamples)*10^2.5;

noise_s = .2*10^4.75;
sl_s    = 5*10^4.75;
bb_s    = 1*10^4.75;


% for generation of pink noise
f       = -n/2:n/2-1;
alpha   = 1;
cpfov   = ifftshift(f);
amp     = 1./(abs(cpfov).^alpha);
amp(1)  = 0;


for ii = 1:num_runs
    fprintf('.');
    for jj = 1:num_stims
        
        bb    = real(ifft(fft(randn(n,num_chan),[],1) .* repmat(amp',[1 num_chan]),[],1));
        noise = real(ifft(fft(randn(n,num_chan),[],1) .* repmat(amp',[1 num_chan]),[],1));
        
        idx = (0:1299)+ triggers((ii-1)*140+jj,1);
        sl      = pred.gain.M_02(jj,:)' * (stim_locked_signal .* sl_s);
        bb      = bsxfun(@times, bb, pred.gain.M_02(jj,:))' *  bb_s;
        noise   = noise' * noise_s;
        
        data157(:,idx) = sl + bb + noise;
    end
end

fprintf('\n')

% data = [data157 data];



cur_time = datestr(now);
cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';

fname = mprfExportDataPath('synthetic',['run_' cur_time '.sqd']);


% sqdwrite(raw_meg_file,fname,data);

fprintf('Copying data file...\n')
system(sprintf('cp %s %s',raw_meg_file,fname))

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
    sqdwrite(fname, fname ,...
        'Channels', data_chan(cur_idx), ...
        'Data',data157(cur_idx,:)',...
        'Action','Overwrite',...
        'Samples',[1 size(data157,2)]);
end

%ft_write_data(fname, data, 'header', hdr, 'dataformat', 'edf');

fname_02 = mprfExportDataPath('synthetic',['timing_info_' cur_time '.mat']);
save(fname_02,'timing','triggers');

cd(cur_dir)

end



