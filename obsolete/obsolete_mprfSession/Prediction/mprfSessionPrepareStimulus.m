function fpath = mprfSessionPrepareStimulus(stim)

if iscell(stim)
    
    load(stim{1});
    load(stim{2});
        
    x_range = [min(grid.Xd(:)) max(grid.Xd(:))];
    y_range = [min(grid.Yd(:)) max(grid.Yd(:))];
    
    stim_range = double([min(stimulus.images(:)) max(stimulus.images(:))]);
    tmp = stimulus.images(:,:,1);
    bk = double(mode(tmp(:)));
        
    new_period = stimulus.seq(stimulus.trigSeq > 0);
    
    im_out_size = [101 101];
    pred_stim.full_im = zeros([im_out_size, length(new_period)]);
    
    srate = max(grid.Xd(:)) ./ floor(max(im_out_size) / 2);
    
    for n = 1:length(new_period)
        tmp_im = double(stimulus.images(:,:,new_period(n)));
        tmp_im = imresize(tmp_im,im_out_size,'nearest');
        pred_stim.full_im(:,:,n) = double(ceil(abs(tmp_im - bk) ./ max(abs(stim_range - bk)))).*(srate.^2);
        
    end
    
    [pred_stim.full_x, pred_stim.full_y] = meshgrid(linspace(x_range(1), x_range(2), im_out_size(1)),...
        linspace(y_range(1), y_range(2), im_out_size(2)));
    
    keep_window = sum(pred_stim.full_im,3) > 0;
    
    pred_stim.im = reshape(pred_stim.full_im,[],size(pred_stim.full_im,3));
    pred_stim.im =  pred_stim.im(keep_window(:),:);
    pred_stim.X = pred_stim.full_x(:);
    pred_stim.Y = pred_stim.full_y(:);
    pred_stim.X = pred_stim.X(keep_window);
    pred_stim.Y = pred_stim.Y(keep_window);
    
elseif isstruct(stim)
    
    if ndims(stim.im) == 2
        pred_stim.im = stim.im;
        
    else
        error('Unexpected stimulus dimensions')
    end
    
    if ndims(stim.X) == 2 && size(stim.X,2) == 1
        pred_stim.X = stim.X;
        
    else
        error('Unexpected grid dimensions')
        
    end
    
    if ndims(stim.Y) == 2 && size(stim.Y,2) == 1
        pred_stim.Y = stim.Y;
        
    else
        error('Unexpected grid dimensions')
        
    end
    
    if ndims(stim.full_im) == 3 && (size(stim.full_im,1) == size(stim.full_im,2))
        pred_stim.full_im = stim.full_im;
        
    else
        error('Unexpected grid dimensions')
        
    end
    
    
    if ndims(stim.full_x) == 2 && (size(stim.full_x,1) == size(stim.full_x,2))
        pred_stim.full_x = stim.full_x;
        
    else
        error('Unexpected grid dimensions')
        
    end
    
    
    if ndims(stim.full_y) == 2 && (size(stim.full_y,1) == size(stim.full_y,2))
        pred_stim.full_y = stim.full_y;
        
    else
        error('Unexpected grid dimensions')
        
    end

    
end


cur_time = datestr(now);
cur_time(cur_time ==':' | cur_time == ' ' | cur_time == '-') = '_';

fpath = mprfExportDataPath('pred_stim',['stimulus_' cur_time]);
save(fpath,'pred_stim');


end