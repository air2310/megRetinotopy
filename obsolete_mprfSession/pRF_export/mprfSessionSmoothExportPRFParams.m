function mprfSessionSmoothExportPRFParams
% This function will export the the pRF parameters we need from the model
% you selected during initilization to a .mat file, located in
% 'subject_dir'/pRF_data/data_dir' and to nifti files, located in
% 'subject_dir'/pRF_data/nifti. Additionally, this function will smooth the
% pRF parameters on the cortical surface and attempts to 'normalize' the
% beta values. The function does not take any inputs

%% To do:
% Weigh parameters by VE while smoothing

%%
% We need the mprfSESSION file:
load('mprfSESSION')
global mprfSESSION

% This function relies on several mrVista functions, so navigate to the
% subject's mrVista directory
cd(mprfSESSION.orig.vista_dir);

% Check if mrVista is on your Matlab path, if not, attempt to add it using
% tbUse.
if exist('niftiRead','file');
    
else
    try
        tbUse('vistasoft');
    catch ME
        error('Vistasoft is not on your Matlab path')
        
    end
    
end

% Check if we still can get to the mprfSESSION code directory:
if exist('mprfLoadRMStimulus_02','file')
    
else
    addpath(mprfSESSION.root.root_path);
    if exist('mprfLoadRMStimulus_02','file')
        
    else
        error('Could not add root directory to path')
        
    end
end
    
    
% Check if all the files we need are exist and are in the expected place:
if ~isfield(mprfSESSION,'orig') || ~isfield(mprfSESSION.orig,'vista_dir')
    error('Cannot find mrVista directory path, mprfSession is not initialized properly')
end

if ~exist(mprfSESSION.orig.vista_dir,'dir')
    error('The vista directory: %s does not exist',mprfSESSION.orig.vista_dir)
end

if ~exist(mprfSESSION.source.ret_model,'file')
    error('No retinotopic model imported')
end

if ~exist(mprfSESSION.source.mrv_class,'file')
    error('No mrVista class file imported')
end

% We need the dataTYPES variable
mrGlobals;

% We need a volume view:
hvol = initHiddenGray;

% Try to work out with data type the retinotopic model was based on by
% matching the directories of the pRF model path to the dataTYPES:
dt = {dataTYPES.name};
dt_match = false(size(dt));

rm_parse = strsplit(mprfSESSION.orig.ret_model,'/');


for n = 1:length(dt);
    cur_dt = dt{n};
    
    if any(strcmpi(cur_dt,rm_parse));
        dt_match(n) = true;
    end
end

% If only one match has been found, we can be confident that we found the
% correct dataTYPE, otherwise, bring up a dialog for the user to select it.
if sum(dt_match) == 1
    use_dt = dt{dt_match};
else
    use_dt = listdlg('listdlg',dt,...
        'SelectionMode','single',...
        'PrompString','Please select the data type to which the retinotopic model belongs');
end

% Close gracefully when no valid data type is selected
if isempty(use_dt)
    return
    
end

% Set the volume view to the current data type and add the RM model
hvol = viewSet(hvol,'curdt',lower(use_dt));
hvol = rmSelect(hvol,1,mprfSESSION.source.ret_model);

% Get the stimulus that was used to fit the current pRF model. This is also
% used later on, so we import this as well
if ~isfield(mprfSESSION.source,'stimulus') || ~exist(mprfSESSION.source.stimulus,'file')
    stim = mprfLoadRMStimulus_02(hvol);
    
    mprfSessionImportSourceData({'rm_stimulus',stim})
    
else
    
    
    
end

% Mask to exclude unreliable voxels (i.e. VE == 0) from the smoothing
% below. Otherwise, the pRF parameters will be averaged with a lot of
% zeros:
sm_mask = rmGet(hvol.rm.retinotopyModels{1},'varexplained') >0;

% We need these parameters from the pRF model
params = {'sigma','x','y','varexplained','beta'};

% We need the mrVista segmentation to check if the selection of pRF
% parameters is correct, i.e. all selected pRF parameters must fall in the
% gray matter.
cls = niftiRead(mprfSESSION.source.mrv_class);

% Initialize the weighted connectivity matrix used in the smoothing
wConMat = [];

% Loop over the needed parameters:
for nn = 1:length(params)
    
    cur_param = params{nn};
    
    % Use mprfExportDataPath to generate a path and file name for the
    % current parameter's NIFTI file:
    fname = mprfExportDataPath('prf_nifti',[cur_param '.nii.gz']);
    
    % Load the current parameter in the VOLUME view. This is mainly used to
    % export the pRF parameters to nifti files. The nifti files are really
    % just for inspection and are not used in any further analyses
    hvol = rmLoad(hvol,1,cur_param,'map');
    hvol = refreshScreen(hvol);
    
    
    % Get the data directly from the retinotopic model as well. This is the
    % data that we use to smooth and store for further analyses:
    if strcmpi(cur_param,'beta')
        tmp = rmGet(hvol.rm.retinotopyModels{1},'bcomp1');
        prf_par_exp.(cur_param) = squeeze(tmp);
        
    else
        
        prf_par_exp.(cur_param) = rmGet(hvol.rm.retinotopyModels{1},cur_param);
    end
    
    % Store the data as nifti and check against the segmentation to see of
    % parameter nifti aligns with the segmentation:
    functionals2nifti(hvol,1 , fname);
    mprfCheckParameterNiftiAlignment(cls, fname);
    
    switch lower(cur_param)
        
        case 'beta'
            
            % Compute the maximum response for every included pRF, by
            % reconstructing the pRF, multiplying it with the stimulus and
            % it's beta, and taking the maximum response from the
            % predicted time series
            mresp = mprfComputeMaximumResponse(stim,sigma_us,x0,y0,prf_par_exp.(cur_param),sm_mask);
            
            % Store the maximum responses as a nifti file:
            fname = mprfExportDataPath('prf_nifti','mresp.nii.gz');
            prf_par_exp.('mresp') = mresp;
            hvol = viewSet(hvol,'map',{mresp});
            functionals2nifti(hvol, 1, fname);
            mprfCheckParameterNiftiAlignment(cls, fname);
            
            prf_par_exp.('mresp') = mresp;
            
            
            % Smooth the maximum responses on the cortical surface
            [mresp_sm, wConMat] = dhkGraySmooth(hvol,mresp,[ ],wConMat, sm_mask);
            
            % Export smoothed maximum responses as a nifti:
            hvol = viewSet(hvol,'map',{mresp_sm});
            fname = mprfExportDataPath('prf_nifti','mresp_smoothed.nii.gz');
            functionals2nifti(hvol, 1, fname);
            mprfCheckParameterNiftiAlignment(cls, fname);
            
            prf_par_exp.('mresp_smoothed') = mresp_sm;
            
            % Recompute the beta using by dividing the smoothed maxumimum
            % responses by the maximum response given the stimulus and
            % smoothed pRF parameters:
            recomp_beta = mprfRecomputeBetas(stim,sigma_smooth,x0_smooth,y0_smooth,mresp_sm);
            
            % Store as nifti:
            fname = mprfExportDataPath('prf_nifti','recomp_beta.nii.gz');
            hvol = viewSet(hvol,'map',{recomp_beta});
            functionals2nifti(hvol, 1, fname);
            mprfCheckParameterNiftiAlignment(cls, fname);
            
            prf_par_exp.('recomp_beta') = recomp_beta;
            
        case 'varexplained'
            
        case {'x','y','sigma'}
            
            % Smooth the current paramter:
            [tmp_sm_par, wConMat] = dhkGraySmooth(hvol,prf_par_exp.(cur_param),[ ],wConMat, sm_mask);
            
            
            % Export smoothed data as nifti:
            fname = mprfExportDataPath('prf_nifti',[cur_param '_smoothed.nii.gz']);
            hvol = viewSet(hvol,'map',{tmp_sm_par});
            functionals2nifti(hvol,1 , fname);
            mprfCheckParameterNiftiAlignment(cls, fname);
            
            prf_par_exp.([cur_param '_smoothed']) = tmp_sm_par;
            
            % Store the current parameters as we need them for the Beta
            % computations above:
            if strcmpi(cur_param,'x')
                x0 = prf_par_exp.(cur_param);
                x0_smooth = tmp_sm_par;
            elseif strcmpi(cur_param,'y')
                y0  = prf_par_exp.(cur_param);
                y0_smooth = tmp_sm_par;
            elseif strcmpi(cur_param,'sigma')
                sigma_us = prf_par_exp.(cur_param);
                sigma_smooth = tmp_sm_par;
            end
            
            
    end
    
end

% WORK ON THIS MASK A BIT MORE.
% It is primarily used for rendering the data on a mesh, but will also
% exclude pRF that are 'obtained' by the smoothing.
prf_par_exp.mask = prf_par_exp.recomp_beta == 0;
fname = mprfExportDataPath('prf_data_file','exported_prf_params.mat');
save(fname, 'prf_par_exp');


cd(mprfSESSION.init.main_dir);
save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
mrvCleanWorkspace;
clear mprfSESSION


end




















