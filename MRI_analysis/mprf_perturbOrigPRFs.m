function mprf_perturbOrigPRFs(prfSurfPath, opt)
% Wrapper function to perturb original pRF parameters on the cortical
% surface.
%       mprf_perturbOrigPRFs(prfSurfPath, dirPth, opt)
% 
% INPUTS:
%   prfSurfPath     : path to surface files containing prf parameters (string)
%   opt             : struct with boolean flags. Should contain the field 
%                     'perturbOrigPRFs' with one of the following definitions:
%                     'position', 'size', 'scramble', or false to exit
%
% OUTPUTS:
%   none
%
%
%
% Author: Eline R. Kupers <ek99@nyu.edu>, 2019
                      
% Select appropiate pertubation function
switch opt.perturbOrigPRFs
    
    case false
        warning('(%s): pRF pertubation option was not defined!', mfilename) 
        return;
    
    case 'position'
        mprf_varyPRFPositionOnSurface(prfSurfPath, opt);
        
    case 'size'
        mprf_varyPRFSizeOnSurface(prfSurfPath, opt);
        
    case 'scramble'
        mprf_scramblePRFOnSurface(prfSurfPath, opt);
end