function [flt_cos, flt_dog] = mprfComputeAndStoreBPFilter(load_file, cpd, bandwidth, fltsz,screenres,imdiam, save_file, do_plot)


if notDefined('bandwidth'); bandwidth = 2; end;
if notDefined('fltsz'); fltsz = 31; end;
if notDefined('screenres'); screenres = 768; end;
if notDefined('cpd'); cpd       = 3; end;
if notDefined('imdiam'); imdiam    = 24; end;
if notDefined('do_plot'); do_plot = false; end;
if notDefined('save_file'); save_file = 'Default'; end;


save_path = '/Volumes/server/Projects/MEG/Retinotopy/Stimuli';
if isempty(load_file)
    
    
    cpim      = cpd * imdiam;
    
    flt_cos = mkBandpassCosine(screenres, cpim, bandwidth, fltsz, 0);
    flt_dog = mkBandpassDog(screenres, cpim, bandwidth, fltsz, 0);
    
    save(fullfile(save_path, 'Filters',save_file),'flt_cos','flt_dog','fltsz','imdiam','screenres','cpd','bandwidth');
    
else
    tmp = load(fullfile(save_path, 'Filters',load_file));
    
    fprintf('Current filter assumes: \n Image diameter = %0.2f degrees\n Screen height = %d pixels\n and is %d pixels wide, centered on %0.2f CPD with %0.2f octaves bandwidth\n',...
        tmp.imdiam, tmp.screenres, tmp.fltsz, tmp.cpd, tmp.bandwidth);
    
    
    if nargin > 1
        if tmp.fltsz == fltsz && ...
                tmp.imdiam == imdiam && ...
                tmp.screenres == screenres && ...
                tmp.cpd == cpd && ...
                tmp.bandwidth == bandwidth
            
            fprintf('Current filter assumes: \n Image diameter = %0.2f degrees\n Screen height = %d pixels\n and is %d pixels wide, centered on %0.2f CPD with %0.2f octaves bandwidth\n',...
                tmp.imdiam, tmp.screenres, tmp.fltsz, tmp.cpd, tmp.bandwidth);
            
            flt_cos = tmp.flt_cos;
            flt_dog = tmp.flt_dog;
            
        else
            fprintf('Requested paramaters do not match the stored parameters, creating new filters');
            
            
            flt_cos = mkBandpassCosine(screenres, cpim, bandwidth, fltsz, 0);
            flt_dog = mkBandpassDog(screenres, cpim, bandwidth, fltsz, 0);
            
            
            
            
        end
    else
        
        flt_cos = tmp.flt_cos;
        flt_dog = tmp.flt_dog;
        
        
        
    end
    
end

if do_plot
    flt = flt_cos;
    
    degperpix = tmp.imdiam/tmp.screenres;
    pixperdeg = tmp.screenres/tmp.imdiam;
    
    x = linspace(0,tmp.imdiam, tmp.screenres+1); x = x(2:end);
    fs = (0:length(x)-1)/max(x);
    
    fltsz_indeg = tmp.fltsz * degperpix;
    x = linspace(0,fltsz_indeg, tmp.fltsz+1);
    x = x(2:end);
    fs = (0:length(x)-1)/max(x);
    
    F = abs(fft2(flt));
    figure, imagesc(fs, fs, F)
    figure, plot(fs, F(:,1))
    
    
    flt = flt_dog;
    
    degperpix = tmp.imdiam/tmp.screenres;
    pixperdeg = tmp.screenres/tmp.imdiam;
    
    x = linspace(0,tmp.imdiam, tmp.screenres+1); x = x(2:end);
    fs = (0:length(x)-1)/max(x);
    
    fltsz_indeg = tmp.fltsz * degperpix;
    x = linspace(0,fltsz_indeg, tmp.fltsz+1);
    x = x(2:end);
    fs = (0:length(x)-1)/max(x);
    
    F = abs(fft2(flt));
    figure, imagesc(fs, fs, F)
    figure, plot(fs, F(:,1))
    
    
end

end