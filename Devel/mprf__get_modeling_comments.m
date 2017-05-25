function comment = mprf__get_modeling_comments(w_model)

comment = 'None';

switch lower(w_model)
    
    case 'run original model'
        comment = 'Generate predictions based on the estimated pRF parameters, using the selected prediction stimulus';
    
    case 'equal weight'
        comment = 'Generate predictions based on the estimated pRF parameters, using the selected prediction stimulus but set the Beta weight to the same value for every pRF';
    
    case 'roi specific predictions'
        comment = 'Generate predictions for each imported ROI specifically';
    
    case 'regress roi predictions'
        comment = 'Regress ROI specific prediction on the measured MEG time series';
    
    case 'scramble prf parameters'
        comment = 'Scramble the pRF parameters between vertices';
        
    case 'fix prf size'
        comment = 'Fix the pRF size to a predefined size';
        
        
    case 'prf size range'
        comment = 'Vary pRF size over a fixed range of pRF sizes';
        
    case 'reliability check'
        
    case 'roi specific'
        
    
        
        
end

















