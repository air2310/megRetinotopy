function out = mprf__compute_meg_measure(data, opts)
out.measure = [];
out.av = [];
out.std = [];
out.ste =[];

% INPUTS:
% current data (i.e. best if we keep this as simple as possible, no loops
% or something. Assumin; ntime_points as first dimension
% What measure to extract
% Do average, if so, return average across second dimension
% Do std, if so, also return standard deviation across second dimension
% Do ste, if so, also return standard error across second dimension
% Just input a options structure that holds these option. Also the

if ~exist('data','var') || isempty(data)
    error('Need data')
end

if ~exist('opts','var') || isempty(opts)
    error('Need options variable for at least the measure type')
end

if ~isfield(opts,'do_av')
    opts.do_av = false;
end


if ~isfield(opts,'do_std')
    opts.do_std = false;
end

if ~isfield(opts,'do_ste')
    opts.do_ste = false;
end

if ~isfield(opts,'do_ft')
    opts.do_ft = true;
end

switch lower(opts.measure)
    
    case {'stim_amp','stim_co','stim_ph','bb_amp','bb_co','bb_ph'}
        
        if ~isfield(opts,'stim_idx')
            error('Need stimulus frequency index')
        end
        
        if opts.do_ft
            ft_data = fft(data);
            ft_data = ft_data(1:1+fix(size(ft_data,1)/2),:);
            sc_amp = abs(ft_data);
        else
            ft_data = data.ft;
            sc_amp = data.sc_amp;
            
            
            
        end
        
        
        if strcmpi(opts.measure,'stim_co') || strcmpi(opts.measure,'bb_co')
            if ~isfield(opts,'inc_idx')
                error('Need indices of remaining frequencies to include')
            end
            sqrtsummagsq = sqrt(sum(sc_amp(opts.inc_idx,:).^2));
            
            
        end
        
        if strcmpi(opts.measure,'stim_ph') || strcmpi(opts.measure,'bb_ph')
            ang_av = @(th) angle(sum(exp(th(:)*1i)));
            ang_std = @(th) sqrt(-2*log((abs(sum(exp(th(:).*1i))) ./ numel(th))));
            
        end
    otherwise
        
        
end

switch lower(opts.measure)
    
    
    case 'stim_amp'
        
        if opts.do_ft
            out.raw = 2*(sc_amp(opts.stim_idx,:))/size(data,1);
        else
            out.raw = 2*(sc_amp(opts.stim_idx,:))/data.size_01;
            
        end
        
        if opts.do_av
            out.av = nanmean(out.raw);
        end
        
        if opts.do_std
            out.std = nanstd(out.raw);
        end
        
        if opts.do_ste
            tmp = nanstd(out.raw);
            out.ste = tmp ./ sqrt(sum(~isnan(tmp)));
        end
        
        
    case 'stim_co'
        out.raw = sc_amp(opts.stim_idx,:) ./ sqrtsummagsq;
        
        
        if opts.do_av
            out.av = nanmean(out.raw);
            
        end
        
        if opts.do_std
            out.std = nanstd(out.raw);
            
        end
        
        if opts.do_ste
            tmp = nanstd(out.raw);
            out.ste = tmp ./ sqrt(sum(~isnan(tmp)));
            
        end
        
    case 'stim_ph'
        
        out.raw = angle(ft_data(sl_idx,:));
        is_nan = isnan(out.raw);
        if opts.do_av
            out.av = ang_av(out.raw(~is_nan)); % Use average angle here, i.e. circular....
            
        end
        
        if opts.do_std
            out.std = ang_std(out.raw(~is_nan));
            
        end
        
        if opts.do_ste
            tmp = ang_std(cur_angle(~is_nan));
            out.ste = tmp./ sqrt(sum(~is_nan));
            
        end
        
        clear is_nan
        
    case 'bb_amp'
        
        
    case 'bb_co'
        
        
    case 'bb_ph'
        
        
        
end



end





























