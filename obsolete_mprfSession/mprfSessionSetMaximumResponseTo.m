function mprfSessionSetMaximumResponseTo(mresp)

if ~exist('cellfind','file')
    tbUse('vistasoft')
end

load(fullfile(pwd,'mprfSESSION.mat'))

load(mprfSESSION.source.rm_stim)

needed_par = {'sigma_smoothed','x_smoothed','y_smoothed'};

load(mprfSESSION.source.bs_head_model);
[~, surf.name] = fileparts(bs_model.SurfaceFile);
surf.type = 'brainstorm';

fname_parts = strsplit(surf.name,'_');

tmp = {'pial','white','mid'};
file_prefix = tmp{ceil(cellfind([strfind(fname_parts,tmp{1}),strfind(fname_parts,tmp{2}),strfind(fname_parts,tmp{3}) ])./length(fname_parts))};
file_prefix = [file_prefix '.'];

for n = 1:length(needed_par)
    
    cur_fname = [file_prefix needed_par{n}];
    cur_path = fullfile(mprfSESSION.prf_exp.bs_surface_data, cur_fname);
    
    cur_data = read_curv(cur_path);
    
    if strcmpi(needed_par{n},'x_smoothed') || strcmpi(needed_par{n},'x')
        X0 = cur_data;
        
    elseif strcmpi(needed_par{n},'y_smoothed') || strcmpi(needed_par{n},'y')
        Y0 = cur_data;
        
        
    elseif strcmpi(needed_par{n},'sigma_smoothed') || strcmpi(needed_par{n},'sigma')
        sigma = cur_data;
        
    else
        error('Unknown data type');
        
        
    end
    

    
end

if numel(mresp) == 1
    mresp = ones(size(sigma)) .* mresp;
end

cur_time = datestr(now);
cur_time(cur_time ==':' | cur_time == ' ' | cur_time == '-') = '_';

recomp_beta = mprfRecomputeBetas(rm_stim,sigma',X0',Y0',mresp');
fname = mprfExportDataPath('prf_bs_surface',[file_prefix 'recomp_beta_' cur_time]);
write_curv(fname,recomp_beta,1);



end


