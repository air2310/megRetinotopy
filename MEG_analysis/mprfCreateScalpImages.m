function out_im = mprfCreateScalpImages(data)


fnames = fieldnames(data);


for nn = 1:length(fnames)
    snr = data.(fnames{nn});
    
    if length(size(snr)) == 3
        snr = nanmean(snr,3);
    end
    
    fh = figure;
    for n = 1:size(snr,1)
        
        cur_data = snr(n,:);
        if all(isnan(cur_data))
            cur_data = ones(size(cur_data));
        end
        
        fh = megPlotMap(cur_data,[0.75 max(snr(:)./0.05) .* 0.05],fh,jet(256),sprintf('Stimulus position %d',n));
        drawnow;
        
        tmp  = getframe(fh);
        
        if n == 1
            out_im.(fnames{nn}) = zeros([size(tmp.cdata) size(snr,1)],'uint8');
            
        end
        out_im.(fnames{nn})(:,:,:,n) = tmp.cdata;
        
        
    end
    
end

end


