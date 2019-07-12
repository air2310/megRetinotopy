function prfs_perturbed = mprf_perturbOrigPRFs(prfSurfPath, dirPth, opt)
% Wrapper function to perturb original pRF parameters on the cortical
% surface.
%       mprf_perturbOrigPRFs(prfSurfPath, dirPth, opt)
% 
% INPUTS:
%   prfSurfPath     : path to surface files containing prf parameters (string)
%   dirPth          : paths to files for given subject (struct)
%   opt             : struct with boolean flags. Should contain the field 
%                     'perturbOrigPRFs' with one of the following definitions:
%                     'position', 'size', 'scramble', or false to exit
%
% OUTPUTS:
%   none
                      
% Select appropiate pertubation function
switch opt.perturbOrigPRFs
    
    case false
        warning('(%s): pRF pertubation option was not defined!', mfilename) 
        return;
    
    case 'position'
        prfs_perturbed = mprf_varyPRFPositionOnSurface(prfSurfPath, opt);
        
    case 'size'
        prfs_perturbed = mprf_varyPRFSizeOnSurface(prfSurfPath, opt);
        
    case 'scramble'
        prfs_perturbed = mprf_scramblePRFOnSurface(prfSurfPath, opt);
end