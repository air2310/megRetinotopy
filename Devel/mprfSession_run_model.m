% % %% TO DO: 
% % - Add options to the GUI that allow to combine ROIs
% dorsal and ventral, left and right 
% - Set the broad band range in an
% intelligent way -> look at Eline's project for defaults 



%%


global mprfSESSION




if exist(fullfile(pwd,'mprfSESSION.mat'),'file')
    load('mprfSESSION.mat');

else
    error('Could not find mprfSESSION file. Please run from session folder')
end


[model, stimulus, syn_data] = mprfSession_model_gui;
[prf, bs, roi] = mprf__model_get_prf_params(model);





return
switch lower(model.type)
    
    
    case 'run original model'
        
        
    case 'equal weight'
        
        
    case 'scramble prf parameters' % NEED SCRAMBLE ITERATIONS
        
        
    case 'fix prf size'
        
        
    case 'pr size range'
        
        
    case 'reliability check'
        
        
        
    case 'fit separate roi predictions'
        
        
        
    otherwise
        
        error('Model type not recognized')
        
        
end






