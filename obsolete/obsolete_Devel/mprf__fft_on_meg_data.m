function ft_data = mprf__fft_on_meg_data(data)

sz = size(data);
n_time = sz(1);
n_bars = sz(2);
n_resp = sz(3);
n_chan = sz(4);

fprintf('Transforming data to frequency domain...')
% Preallocate output:
ft_data = nan(1+fix(n_time/2), n_bars, n_resp, n_chan);

% Loop over channels
for this_chnl = 1:n_chan
    % Loop over stimulus frames
    for this_bar = 1:n_bars
        % FFT
        tmp = fft(data(:,this_bar,:,this_chnl));
        % Remove redundant part
        ft_data(:,this_bar,:,this_chnl) = tmp(1:1+fix(n_time/2),:,:,:);
    end
end

fprintf('Done.\n')

end

