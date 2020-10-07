function out_path = mprfExportDataPath(type, fname)

global mprfSESSION
add_fname = true;



switch lower(type)
    
    case 'prf_nifti'
        
        if isfield(mprfSESSION,'prf_exp') && isfield(mprfSESSION.prf_exp,'nifti_dir')
            
        else
            
            if ~exist(fullfile(mprfSESSION.init.main_dir,'prf_data'),'dir')
                mkdir(mprfSESSION.init.main_dir,'prf_data')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data'),'nifti')
            elseif ~exist(fullfile(mprfSESSION.init.main_dir,'prf_data','nifti'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data'),'nifti')
            end
            mprfSESSION.prf_exp.nifti_dir = fullfile(mprfSESSION.init.main_dir,'prf_data','nifti');
            save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
        end
        
        out_path = mprfSESSION.prf_exp.nifti_dir;
        
    case 'prf_data_file'
        if isfield(mprfSESSION,'prf_exp') && isfield(mprfSESSION.prf_exp,'data_dir')
            
        else
            
            if ~exist(fullfile(mprfSESSION.init.main_dir,'prf_data'),'dir')
                mkdir(mprfSESSION.init.main_dir,'prf_data')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data'),'data_dir')
            elseif ~exist(fullfile(mprfSESSION.init.main_dir,'prf_data','data_dir'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data'),'data_dir')
            end
            mprfSESSION.prf_exp.data_dir = fullfile(mprfSESSION.init.main_dir,'prf_data','data_dir');
            save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
            
        end
        
        out_path = mprfSESSION.prf_exp.data_dir;
        
        
        
    case {'prf_fs_surface','wang_prf_fs_surface'}
        if isfield(mprfSESSION,'prf_exp') && isfield(mprfSESSION.prf_exp,'fs_surface_data')
            
        else
            if ~exist(fullfile(mprfSESSION.init.main_dir,'prf_data'),'dir')
                mkdir(mprfSESSION.init.main_dir,'prf_data')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data'),'surface')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data','surface'),'freesurfer')
                
            elseif ~exist(fullfile(mprfSESSION.init.main_dir,'prf_data','surface'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data'),'surface')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data','surface'),'freesurfer')
                
            elseif ~exist(fullfile(mprfSESSION.init.main_dir,'prf_data','surface','freesurfer'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data','surface'),'freesurfer')
            end
            mprfSESSION.prf_exp.fs_surface_data = fullfile(mprfSESSION.init.main_dir,'prf_data','surface','freesurfer');
            save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
            
        end
        
        out_path = mprfSESSION.prf_exp.fs_surface_data;
        
        
        
        
    case 'prf_bs_surface'
        if isfield(mprfSESSION,'prf_exp') && isfield(mprfSESSION.prf_exp,'bs_surface_data')
            
        else
            if ~exist(fullfile(mprfSESSION.init.main_dir,'prf_data'),'dir')
                mkdir(mprfSESSION.init.main_dir,'prf_data')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data'),'surface')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data','surface'),'brainstorm')
                
            elseif ~exist(fullfile(mprfSESSION.init.main_dir,'prf_data','surface'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data'),'surface')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data','surface'),'brainstorm')
                
            elseif ~exist(fullfile(mprfSESSION.init.main_dir,'prf_data','surface','brainstorm'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'prf_data','surface'),'brainstorm')
            end
            mprfSESSION.prf_exp.bs_surface_data = fullfile(mprfSESSION.init.main_dir,'prf_data','surface','brainstorm');
            save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
            
        end
        
        out_path = mprfSESSION.prf_exp.bs_surface_data;
        
        
    case {'roi_fs_surface','wang_roi_fs_surface'}
        
        if isfield(mprfSESSION,'roi_exp') && isfield(mprfSESSION.roi_exp,'freesurfer')
            
        else
            if ~exist(fullfile(mprfSESSION.init.main_dir,'rois'),'dir')
                mkdir(mprfSESSION.init.main_dir,'rois')
                mkdir(fullfile(mprfSESSION.init.main_dir,'rois'),'surface')
                mkdir(fullfile(mprfSESSION.init.main_dir,'rois','surface'),'freesurfer')
                
            elseif ~exist(fullfile(mprfSESSION.init.main_dir,'rois','surface'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'rois'),'surface')
                mkdir(fullfile(mprfSESSION.init.main_dir,'rois','surface'),'freesurfer')
                
            elseif ~exist(fullfile(mprfSESSION.init.main_dir,'rois','surface','freesurfer'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'rois','surface'),'freesurfer')
            end
            mprfSESSION.roi_exp.freesurfer = fullfile(mprfSESSION.init.main_dir,'rois','surface','freesurfer');
            save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
            
        end
        
        out_path = mprfSESSION.roi_exp.freesurfer;
        
        
        
    case 'roi_bs_surface'
        
        if isfield(mprfSESSION,'roi_exp') && isfield(mprfSESSION.roi_exp,'brainstorm')
            
        else
            if ~exist(fullfile(mprfSESSION.init.main_dir,'rois'),'dir')
                mkdir(mprfSESSION.init.main_dir,'rois')
                mkdir(fullfile(mprfSESSION.init.main_dir,'rois'),'surface')
                mkdir(fullfile(mprfSESSION.init.main_dir,'rois','surface'),'brainstorm')
                
            elseif ~exist(fullfile(mprfSESSION.init.main_dir,'rois','surface'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'rois'),'surface')
                mkdir(fullfile(mprfSESSION.init.main_dir,'rois','surface'),'brainstorm')
                
            elseif ~exist(fullfile(mprfSESSION.init.main_dir,'rois','surface','brainstorm'),'dir')
                mkdir(fullfile(mprfSESSION.init.main_dir,'rois','surface'),'brainstorm')
            end
            mprfSESSION.roi_exp.brainstorm = fullfile(mprfSESSION.init.main_dir,'rois','surface','brainstorm');
            save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
            
        end
        
        out_path = mprfSESSION.roi_exp.brainstorm;
        
    case 'pred_stim'
        
        if ~exist(fullfile(mprfSESSION.init.main_dir,'prediction'),'dir');
            mkdir(mprfSESSION.init.main_dir,'prediction');
            
        end
        
        if ~exist(fullfile(mprfSESSION.init.main_dir,'prediction','stimuli'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'prediction'),'stimuli');
            
        end
        
        
        if ~isfield(mprfSESSION,'pred') || ~isfield(mprfSESSION.pred,'stimuli')
            mprfSESSION.pred.stimuli = cell(1,1);
            mprfSESSION.pred.stimuli{1,1} = fullfile(mprfSESSION.init.main_dir,'prediction','stimuli',fname);
            
        else
            mprfSESSION.pred.stimuli{end+1} = fullfile(mprfSESSION.init.main_dir,'prediction','stimuli',fname);
            
        end
        
        save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
        
        out_path = mprfSESSION.pred.stimuli{end};
        
        add_fname = false;
        
        
    case 'prediction'
        if ~exist(fullfile(mprfSESSION.init.main_dir,'prediction'),'dir');
            mkdir(mprfSESSION.init.main_dir,'prediction');
            
        end
        
        if ~exist(fullfile(mprfSESSION.init.main_dir,'prediction','MEG_prediction_runs'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'prediction'),'MEG_prediction_runs');
            
        end
        
        if ~isfield(mprfSESSION,'pred') || ~isfield(mprfSESSION.pred,'runs')
            mprfSESSION.pred.runs = cell(1,1);
            mprfSESSION.pred.runs{1,1} = fullfile(mprfSESSION.init.main_dir,'prediction','MEG_prediction_runs',fname);
            
        else
            mprfSESSION.pred.stimuli{end+1} = fullfile(mprfSESSION.init.main_dir,'prediction','MEG_prediction_runs',fname);
            
        end
        
        save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
        out_path = mprfSESSION.pred.runs{end};
        
        add_fname = false;
        
    case 'synthetic'
         if ~exist(fullfile(mprfSESSION.init.main_dir,'prediction'),'dir');
            mkdir(mprfSESSION.init.main_dir,'prediction');
            
        end
        
        if ~exist(fullfile(mprfSESSION.init.main_dir,'prediction','synthetic'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'prediction'),'synthetic');
            
        end
        
        if ~isfield(mprfSESSION,'pred') || ~isfield(mprfSESSION.pred,'syn')
            mprfSESSION.pred.syn = fullfile(mprfSESSION.init.main_dir,'prediction','synthetic');
        end
        
        save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
        out_path = mprfSESSION.pred.syn;
        
    case 'pred_gain'
        if ~exist(fullfile(mprfSESSION.init.main_dir,'prediction'),'dir');
            mkdir(mprfSESSION.init.main_dir,'prediction');
            
        end
        
        if ~exist(fullfile(mprfSESSION.init.main_dir,'prediction','gain'),'dir')
            mkdir(fullfile(mprfSESSION.init.main_dir,'prediction'),'gain');
            
        end
        
        if ~isfield(mprfSESSION,'pred') || ~isfield(mprfSESSION.pred,'gain')
            mprfSESSION.pred.gain = cell(1,1);
            mprfSESSION.pred.gain{1,1} = fullfile(mprfSESSION.init.main_dir,'prediction','gain',fname);
            
        else
            mprfSESSION.pred.gain{end+1} = fullfile(mprfSESSION.init.main_dir,'prediction','gain',fname);
            
        end
        
        save(fullfile(mprfSESSION.init.main_dir,'mprfSESSION'),'mprfSESSION')
        out_path = mprfSESSION.pred.gain{end};
        
        add_fname = false;
        
        
        
    otherwise
        
        
        error('Unknown type')
        
end

if add_fname
    out_path = fullfile(out_path,fname);
    
else
    
    
end




end































