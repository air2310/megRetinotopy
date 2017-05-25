function full_roi_tag = mprf__get_roi_tags(roi_in)


prim_rnames = {'v1','v2','v3','hv4','lo1','lo2','vo1','vo2','v3ab','to1','to2'};
sec_rnames= {'v4','v3a','v3b'};
dv_tags = {'d','v'};
lr_tags = {'l','r'};

[~, cur_roi] = fileparts(lower(roi_in));

roi_tag = '';
dv_tag = '';
lr_tag = '';
ask_for_tag = true;

prim_match = [];
for nn = 1:length(prim_rnames)
    cur_lbl = prim_rnames{nn};
    
    found = [regexp(cur_roi,['.' cur_lbl '.'])...
        regexp(cur_roi,['\<' cur_lbl]) ...
        regexp(cur_roi,[cur_lbl '\>' ])];
    if ~isempty(found) && length(found) == 1
        prim_match = [prim_match nn];
        found = [];
        
        
    end
    
end



sec_match = [];
for nn = 1:length(sec_rnames)
    cur_lbl = sec_rnames{nn};
    
    found = [regexp(cur_roi,['.' cur_lbl '.'])...
        regexp(cur_roi,['\<' cur_lbl]) ...
        regexp(cur_roi,[cur_lbl '\>' ])];
    if ~isempty(found) && length(found) == 1
        sec_match = [sec_match nn];
        found = [];
        
        
    end
    
end

if ~isempty(prim_match) && ~isempty(sec_match)
    roi_tag_01 = prim_rnames{prim_match};
    roi_tag_02 = sec_rnames{sec_match};
    
    if length(roi_tag_01) > length(roi_tag_02)
        roi_tag = roi_tag_01;
        
    elseif length(roi_tag_01) < length(roi_tag_02)
        roi_tag = roi_tag_02;
        
    elseif length(roi_tag_01) == length(roi_tag_02)
        
        
        
    end
    
elseif ~isempty(prim_match) && isempty(sec_match)
    if length(prim_match) == 1
        m_idx = strfind(cur_roi,prim_rnames{prim_match});
        cur_roi(m_idx : m_idx+length(prim_rnames{prim_match})-1) = [];
        roi_tag = prim_rnames{prim_match};
    end
    
elseif isempty(prim_match) && ~isempty(sec_match)

    if length(sec_match) == 1
        m_idx = strfind(cur_roi,prim_rnames{prim_match});
        
        cur_roi(m_idx : m_idx+length(prim_rnames{prim_match})-1) = [];
        roi_tag = sec_rnames{sec_match};
    end
    
end



if isempty(roi_tag)
    
    
else
    dv_match = [];
    for nn = 1:length(dv_tags)
        
        found = regexp(cur_roi,[dv_tags{nn} '\>' ]);
        if ~isempty(found) && length(found) == 1
            dv_match = [dv_match nn];
            found = [];
            
        end
    end
    
    if ~isempty(dv_match)
        if dv_match == 1
            dv_tag = 'd';
            
        elseif dv_match == 2
            dv_tag = 'v';
            
            
            
        end
        
    end
    
    lr_match = [];
    for nn = 1:length(lr_tags)
        
        found = regexp(cur_roi,['\<' lr_tags{nn} ]);
        if ~isempty(found) && length(found) == 1
            lr_match = [lr_match nn];
            found = [];
            
        end
        if ~isempty(lr_match)
            if lr_match == 1
                lr_tag = 'left';
                
            elseif lr_match == 2
                lr_tag = 'right';
                
            end
            
        end
        
    end
end


if isempty(lr_tag) || isempty(roi_tag)
    
    
    
else
    if isempty(dv_tag)
        full_roi_tag = [lr_tag '_' roi_tag];
    else
        full_roi_tag = [lr_tag '_' roi_tag '_' dv_tag];
    end
    
    ask_for_tag = false;
end

if ask_for_tag
    full_roi_tag = inputdlg(sprintf('Could not match %s to a tag, please provide tag',roi_in));
end
end