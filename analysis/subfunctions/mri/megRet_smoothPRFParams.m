function megRet_smoothPRFParams(s)

% INPUTS:
% s         : struct with paths to important files (created by
%               s_initSession_wlsubjXXX)

% Navigate to the Vista soft session
mainDir = pwd;
cd(s.vistaSession.pth)

% Load need the dataTYPES variable
mrGlobals;

% We need a volume view:
hvol = initHiddenGray;

% Try to work out with data type the retinotopic model was based on by
% matching the directories of the pRF model path to the dataTYPES:
dt = {dataTYPES.name};

for n = 1:length(dt)
    cur_dt = dt{n};
    
    if any(strcmpi(cur_dt,'Averages'))
        break;
    end
end


% Set the volume view to the current data type and add the RM model
hvol = viewSet(hvol,'curdt',cur_dt);
hvol = rmSelect(hvol,1,s.PRFParams.pth);


% Mask to exclude unreliable voxels (i.e. VE == 0) from the smoothing
% below. Otherwise, the pRF parameters will be averaged with a lot of
% zeros:
VEmask = rmGet(hvol.rm.retinotopyModels{1},'varexplained') >0;

% We need these parameters from the pRF model
paramLabels = {'sigma','x','y','varexplained','beta'};

% We need the mrVista segmentation to check if the selection of pRF
% parameters is correct, i.e. all selected pRF parameters must fall in the
% gray matter.
classFile = niftiRead(fullfile(s.vistaSession.pth,'3DAnatomy','t1_class.nii.gz'));
load(fullfile(s.PRFParams.pth));
% Initialize the weighted connectivity matrix used in the smoothing
wConMat = [];

savePth = fullfile(s.outPut.pth, 'prfParams');
if ~exist('savePth', 'dir'); mkdir(savePth); end

% Loop over the needed parameters:
for nn = 1:length(paramLabels)
    
    % Select parameter label
    curParam = paramLabels{nn};
    
    % Generate file name for the current parameter's NIFTI file:
    fname = fullfile(savePth,[curParam '.nii.gz']);
    
    % Load the current parameter in the VOLUME view. This is mainly used to
    % export the pRF parameters to nifti files. The nifti files are really
    % just for inspection and are not used in any further analyses
    hvol = rmLoad(hvol,1,curParam,'map');
    hvol = refreshScreen(hvol);
    
    % Get the data directly from the retinotopic model as well. This is the
    % data that we use to smooth and store for further analyses:
    if strcmpi(curParam,'beta') % [EK]: why not beta?
        tmp = rmGet(hvol.rm.retinotopyModels{1},'bcomp1');
        prfParamsExp.(curParam) = squeeze(tmp);
        
    else        
        prfParamsExp.(curParam) = rmGet(hvol.rm.retinotopyModels{1},curParam);
    end
    
    % Store the data as nifti and check against the segmentation to see of
    % parameter nifti aligns with the segmentation:
    functionals2nifti(hvol,1 , fname);
    mprfCheckParameterNiftiAlignment(classFile, fname);
    
    switch lower(curParam)
        
        case 'beta'
            % Compute the maximum response for every included pRF, by
            % reconstructing the pRF, multiplying it with the stimulus and
            % it's beta, and taking the maximum response from the
            % predicted time series
            mresp = mprfComputeMaximumResponse(s.stim.MRI,sigma_us,x0,y0,prfParamsExp.(curParam),VEmask);
            
            % Store the maximum responses as a nifti file:
            fname = fullfile(savePth,'mresp.nii.gz');
            prfParamsExp.('mresp') = mresp;
            hvol = viewSet(hvol,'map',{mresp});
            functionals2nifti(hvol, 1, fname);
            mprfCheckParameterNiftiAlignment(classFile, fname);
            
            prfParamsExp.('mresp') = mresp;
            
            
            % Smooth the maximum responses on the cortical surface
            [mresp_sm, wConMat] = dhkGraySmooth(hvol,mresp,[ ],wConMat, VEmask);
            
            % Export smoothed maximum responses as a nifti:
            hvol = viewSet(hvol,'map',{mresp_sm});
            fname = fullfile(savePth, 'mresp_smoothed.nii.gz');
            functionals2nifti(hvol, 1, fname);
            mprfCheckParameterNiftiAlignment(classFile, fname);
            
            prfParamsExp.('mresp_smoothed') = mresp_sm;
            
            % Recompute the beta using by dividing the smoothed maxumimum
            % responses by the maximum response given the stimulus and
            % smoothed pRF parameters:
            recomp_beta = mprfRecomputeBetas(s.stim.MRI,sigma_smooth,x0_smooth,y0_smooth,mresp_sm);
            
            % Store as nifti:
            fname = fullfile(savePth,'recomp_beta.nii.gz');
            hvol = viewSet(hvol,'map',{recomp_beta});
            functionals2nifti(hvol, 1, fname);
            mprfCheckParameterNiftiAlignment(classFile, fname);
            
            prfParamsExp.('recomp_beta') = recomp_beta;
            
        case 'varexplained'
            
        case {'x','y','sigma'}
            
            % Smooth the current paramter:
            [tmp_sm_par, wConMat] = dhkGraySmooth(hvol,prfParamsExp.(curParam),[ ],wConMat, VEmask);
            
            
            % Export smoothed data as nifti:
            fname = fullfile(savePth,[curParam '_smoothed.nii.gz']);
            hvol = viewSet(hvol,'map',{tmp_sm_par});
            functionals2nifti(hvol,1 , fname);
            mprfCheckParameterNiftiAlignment(classFile, fname);
            
            prfParamsExp.([curParam '_smoothed']) = tmp_sm_par;
            
            % Store the current parameters as we need them for the Beta
            % computations above:
            if strcmpi(curParam,'x')
                x0 = prfParamsExp.(curParam);
                x0_smooth = tmp_sm_par;
            elseif strcmpi(curParam,'y')
                y0  = prfParamsExp.(curParam);
                y0_smooth = tmp_sm_par;
            elseif strcmpi(curParam,'sigma')
                sigma_us = prfParamsExp.(curParam);
                sigma_smooth = tmp_sm_par;
            end                 
    end % switch    
end % paramLabels

% It is primarily used for rendering the data on a mesh, but will also
% exclude pRF that are 'obtained' by the smoothing.
prfParamsExp.mask = prfParamsExp.recomp_beta == 0;
fname = fullfile(savePth,'allPRFParams.mat');
save(fname, 'prfParamsExp');

cd(mainDir);
mrvCleanWorkspace;