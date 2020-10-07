function data = mprfFilterDataHighPass(data, Fstop, Fpass, Astop, Apass, srate)


if ~exist('Fstop','var') || isempty(Fstop)
    Fstop = 0.1;    %Hz
end

if ~exist('Fpass','var') || isempty(Fpass)
    Fpass = 1;      % Hz
    
end

if ~exist('Astop','var') || isempty(Astop)
    Astop = 60;     % Db
end

if ~exist('Apass','var') || isempty(Apass)
    Apass = 3;      % Passband Ripple (Db)
end

if ~exist('srate','var') || isempty(srate)
    srate = 1000;   % Sampling rate
end

h = fdesign.highpass(Fstop, Fpass, Astop, Apass, srate);
Hd = design(h,'butter','matchexactly','stopband');

fvtool(Hd);

orig_dim = size(data);

fprintf('Filtering data from %d channels:\n',orig_dim(4))
for n = 1:size(data,4)
    fprintf('%d.',n);

    tmp_data = data(:,:,:,n);
    tmp_data = tmp_data(:);
    nan_idx = isnan(tmp_data);
    
    filtered_data = nan(size(tmp_data));
    filtered_data(~nan_idx) = filter(Hd,tmp_data(~nan_idx));

    data(:,:,:,n) = reshape(filtered_data,[orig_dim(1:3) 1]);

end

fprintf('\nDone.\n');

end


% Fstop1 = 0.1;
% Fstop2 = 60;
% Fpass1 = 1;
% Fpass2 = 40;
% Astop1 = 60;
% Astop2 = 60;
% Apass1 = 3;
% Apass2 = 3;
% 
% 
% h_bp = fdesign.bandpass(Fstop1, Fpass1, Fpass2, Fstop2, Astop1, Apass1, Astop2,srate);
% Hd_bp = design(h_bp, 'butter','matchexactly',match);
% 
% fvtool(Hd_bp)









