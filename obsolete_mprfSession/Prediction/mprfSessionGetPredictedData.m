function [stim_im,head_im,M_02,M]  = mprfSessionGetPredictedData

% This is getting a bit messy, this file creates predictions that are used
% later on (by the create synthetic data for example), but also makes head
% images of the predictions, which takes very long to load. Change this
% file to generate the predictions only, or both.

global mprfSESSION
if isempty(mprfSESSION)
    
    load(fullfile(pwd,'mprfSESSION.mat'));
    global mprfSESSION
end



make_predictions = true;

if  isfield(mprfSESSION.pred,'runs')
    
    fnames = cell(size(mprfSESSION.pred.runs));
    for n = 1:length(fnames)
        [~, fnames{n}] = fileparts(mprfSESSION.pred.runs{n});
        
        
        
    end
    answer = listdlg('ListString',fnames,...
        'SelectionMode','single',...
        'PromptString','Please select predictions to display',...
        'CancelString','Make new prediction');
    
    
    if isempty(answer)
    else
        
        pred_to_load = mprfSESSION.pred.runs{answer};
    end
    try
        load(pred_to_load)
        assert(exist('stim_im','var')==1)
        assert(exist('head_im','var')==1)
        assert(exist('M','var')==1)
        assert(exist('M_02','var')==1)
        make_predictions = false;
        
    catch
        fprintf('Failed to load predictions, creating new ones');
        
    end
    
end

if make_predictions
    load(mprfSESSION.source.bs_head_model);
    
    [~, surf.name] = fileparts(bs_model.SurfaceFile);
    surf.type = 'brainstorm';
    
    [pred_resp, rm_stim] = mprfSessionGeneratePRFResponsePredictions(surf, [],'meg');
    
    G = bst_gain_orient(bs_model.Gain,bs_model.GridOrient);
    G2 = G(~isnan(G(:,1)),:); % Remove NaNs....
    keep_sensors = ~isnan(G(:,1));
    
    M = zeros(size(rm_stim.im,2), size(G,1));
    M(:,keep_sensors) = pred_resp * G2';
    
    
    M_02 = M(:,1:157);
    d_range = [min(M_02(:)) max(M_02(:))];
    
    fh_ft = figure;
    
    for n = 1:size(M_02,1);
        
        fh_ft = megPlotMap(M_02(n,:),d_range,fh_ft,jet(256),sprintf('Stimulus position %d',n));
        drawnow;
        
        tmp  = getframe(fh_ft);
        
        if n == 1
            head_im = zeros([size(tmp.cdata) size(M_02,1)],'uint8');
            
        end
        head_im(:,:,:,n) = tmp.cdata;
        
    end
    
    
    cur_time = datestr(now);
    cur_time(cur_time == ' ' | cur_time == ':' | cur_time == '-') = '_';
    
    stim_im = rm_stim.full_im; %#ok<NASGU>
    fname = mprfExportDataPath('prediction',['run_' cur_time]);
    save(fname,'stim_im','head_im','M_02','M');
    
end




end