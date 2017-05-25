function mprfSessionFillPathVariable

global mprfSESSION
global paths

fnames = fieldnames(mprfSESSION.orig);

for n = fnames'
    
    if strcmp(n{1},'vista_dir')
        paths.orig.vista_dir = mprfSESSION.orig.(n{1});
        
    elseif strcmp(n{1},'ret_model')
        paths.orig.rm_file =  mprfSESSION.orig.(n{1});
        
     elseif strcmp(n{1},'fs_subject_dir')
        paths.orig.fs_subject_dir = mprfSESSION.orig.(n{1});
        
    elseif strcmp(n{1},'mr_vista_t1')
        paths.orig.mrvt1  = mprfSESSION.orig.(n{1});
        
     elseif strcmp(n{1},'mr_vista_class')
        paths.orig.mrv_class  = mprfSESSION.orig.(n{1});
        
    elseif strcmp(n{1},'')
        paths.orig.bs_model_file = mprfSESSION.orig.(n{1});
        
    elseif strcmp(n{1},'')
        paths.orig.bs_anat_dir  = mprfSESSION.orig.(n{1});
        
    end
    
    
    
    
    
end








end








