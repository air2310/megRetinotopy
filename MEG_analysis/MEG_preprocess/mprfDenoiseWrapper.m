function [epoch, results, evalout,denoised_spec] = mprfDenoiseWrapper(epoch, stim_dir)


skip_n_samples = 0;
keep_n_samples = 1100;
fs_range = [10 10];

epoch_idx = epoch.idx;
epoch_idx = epoch_idx(:,ones(1,19));
epoch_idx = epoch_idx(:);

first_idx = epoch.preproc.first_pos_idx;
first_idx = first_idx(:,ones(1,19));
first_idx = first_idx(:);

blink_idx = epoch.preproc.blink_idx;
blink_idx = blink_idx(:,ones(1,19));
blink_idx = blink_idx(:);

epoch_idx = epoch_idx(~first_idx &~blink_idx);
epoch_idx = epoch_idx(~epoch.preproc.badEpochs);

data_preproc = epoch.data;

epoch = rmfield(epoch,'data');

fprintf('Reshaping input variable...\n')
data_preproc = reshape(data_preproc,size(data_preproc,1), [], size(data_preproc,4));
data_preproc = permute(data_preproc, [3 1 2]);

fprintf('Removing %d samples from beginning of each period, keeping %d subsequent samples...\n',...
    skip_n_samples, keep_n_samples);
data_preproc = data_preproc(:,skip_n_samples+1:keep_n_samples+skip_n_samples,:);
fprintf('Done.\n')

fprintf('Removing bad epochs and channels and blink and first bar position periods...\n')
data_preproc = data_preproc(:,:,~first_idx & ~blink_idx);
data_preproc = data_preproc(~epoch.preproc.badChannels,:,~epoch.preproc.badEpochs);
fprintf('Done.\n')

all_cond = 1:8;

presented_conds = all_cond(ismember(all_cond,unique(epoch_idx)));

fprintf('Making design matrix...\n')
fprintf('Loading stimulus...\n')
stimfiles = dir(fullfile(stim_dir,'*.mat'));
stim = load(fullfile(stim_dir, stimfiles(2).name));

stim_range = double([min(stim.stimulus.images(:)) max(stim.stimulus.images(:))]);
bk = double(mode(mode(stim.stimulus.images(:,:,1))));
new_period = stim.stimulus.seq(stim.stimulus.trigSeq > 0);

im_out_size = [101 101];
im_seq = zeros([prod(im_out_size), length(new_period)],'uint8');


for n = 1:length(new_period)
    tmp_im = double(stim.stimulus.images(:,:,new_period(n)));
    tmp_im = imresize(tmp_im,im_out_size,'nearest');
    im_seq(:,n) = uint8(ceil(abs(tmp_im(:) - bk) ./ max(abs(stim_range - bk))));
    
end

design_02 = zeros(size(epoch.idx));
n_cond = 0;
for n = 1:length(design_02)
    
    if any(ismember(presented_conds, epoch.idx(n))) ...
            && design_02(n) == 0;
        n_cond = n_cond + 1;
        
        design_02(n) = n_cond;
        
        
        
        for nn = n+1:length(design_02)
            
            if isequal(im_seq(:,n),im_seq(:,nn)) && design_02(nn) == 0;
                design_02(nn) = n_cond;
                
                
            end
        end
    end
    
    
end

design_02 = design_02(:,ones(1,19));
design_02 = design_02(:);

design_02 = design_02(~first_idx &~blink_idx);
design_02 = design_02(~epoch.preproc.badEpochs);

n_conds = max(design_02);
design_03 = zeros(size(design_02,1),n_conds);
cond_idx = sub2ind(size(design_03),find(design_02),design_02(design_02>0));
design_03(cond_idx) = 1;

fprintf('Desing has %d conditions, corresponding to unique bar positions\n',size(design_03,2));

opt.verbose             = true;
opt.use3Channels        = false;
opt.removeFirstEpoch    = false;
opt.removeMsEpochs      = false;
opt.pcchoose          = 1.05;   % Take 10 PCs
opt.npcs2try          = 10;

evokedfun = @(x)mprfDenoiseEvalFun(x,fs_range,1000);

fprintf('Denoising data...\n')

[results,evalout,denoised_spec, denoised_data] = denoisedata(design_03,data_preproc,evokedfun,evokedfun,opt);

fprintf('Done.\n')

figure;
plot(results.finalmodel.r2,results.origmodel.r2,'k.');
hold on;
plot([0 10],[0 10],'r-');
xlabel('Final');
ylabel('orig');

orig_dim = [size(denoised_data{1},2) 2660 157];

first_pos_idx = epoch.preproc.first_pos_idx;
blink_idx = epoch.preproc.blink_idx;

first_pos_idx = first_pos_idx(:,ones(1,19));
blink_idx = blink_idx(:,ones(1,19));

first_pos_idx = first_pos_idx(:);
blink_idx = blink_idx(:);

all_idx = first_pos_idx | blink_idx;

good_epoch_idx = false(size(all_idx));
good_epoch_idx(~all_idx) = ~epoch.preproc.badEpochs;

fprintf('Preallocating output variable.\n Estimated size: %d bites\n',...
    prod([orig_dim 8]))
tmp_data = nan(orig_dim);

fprintf('Done.\n')

tmp_data(:,good_epoch_idx,~epoch.preproc.badChannels) = permute(denoised_data{1},[2 3 1]);

epoch.data = reshape(tmp_data,[1100 140 19 157]);


end


























