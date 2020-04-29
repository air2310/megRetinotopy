function [epoched_data, fpath] = mprfGetGuiData


if ~exist('raw_data','var') || isempty(raw_data)

    [data_fname, data_fpath] = uigetfile('*','Please select data to load');
    
    tmp = load(fullfile(data_fpath,data_fname));
    
    tmp_name = fieldnames(tmp);
    
    epoched_data.data = tmp.(tmp_name{1}).data;
    epoched_data.start_end = tmp.(tmp_name{1}).start_end;
    epoched_data.idx = tmp.(tmp_name{1}).idx;
    epoched_data.preproc = tmp.(tmp_name{1}).preproc;
    fpath = fullfile(data_fpath, data_fname);
    
    return
end

end
% 
% 
% if ischar(raw_data)
%     load(raw_data);
%     raw_data = data_chan_data;
%     
%     data_dir = fileparts(raw_data);
%     
% end
% 
% if isnumeric(raw_data)
% 
%     if ~exist('epoch_length','var') || isempty(epoch_length)
%         epoch_length = 'full';
%         
%     end
%     
%     
%     if ~exist('triggers','var') || isempty(triggers)
%         error('Need triggers variable');
%     
%         
%     end
%     
%     if ~exist('do_save','var')
%         do_save = [];
%     end
%     
%     
%     if ~exist('data_dir','var')
%         data_dir = '';
%         
%     end
%     
%     
%     epoched_data = mprfEpochMEGData(raw_data, triggers, epoch_length ,do_save, data_dir);
% 
% 
% 
% end
% 
% 
% 
% 
