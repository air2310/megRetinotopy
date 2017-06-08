function bs = mprf__get_lead_field(bs)


if isempty(which('bst_gain_orient.m'))
    error('Need brainstorm files')
end

load(bs.model_file)

bs.lead_field = bst_gain_orient(bs_model.Gain,bs_model.GridOrient);
bs.lead_field2 = bs.lead_field(~isnan(bs.lead_field(:,1)),:); % Remove NaNs....
bs.keep_sensors = ~isnan(bs.lead_field(:,1));


%
% fname_parts = strsplit(bs_model.SurfaceFile,'_');
%
% tmp = {'pial','white','mid'};
%
% file_prefix = tmp{ceil(find(~cellfun(@isempty, [cellfun(@(x) strfind(x,tmp{1}), fname_parts, 'UniformOutput',false) ...
%     cellfun(@(x) strfind(x,tmp{2}), fname_parts, 'UniformOutput',false) ...
%     cellfun(@(x) strfind(x,tmp{3}), fname_parts, 'UniformOutput',false)])) ...
%     ./length(fname_parts))};
%

end




