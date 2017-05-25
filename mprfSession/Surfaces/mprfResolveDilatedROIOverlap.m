function pruned_inds = mprfResolveDilatedROIOverlap(inds_dil, inds_orig)

pruned_inds = cell(size(inds_dil));

for n = 1:length(inds_dil)
    dil_inds = unique(inds_dil{n});
    orig_inds2 = unique(inds_orig{n});
    
    added = ~ismember(dil_inds, orig_inds2);
    added_inds = dil_inds(added);
    
    for nn = 1:length(inds_orig)
        orig_inds = unique(inds_orig{nn});
        
        if n == nn
            % Skip, comparing the same ROIs
            
        else
            overlap = ismember(added_inds,orig_inds);
            added_inds = added_inds(~overlap);
        end
        
        
        
    end
    pruned_inds{n} = unique([added_inds(:); orig_inds2(:)]);
    
end

end










