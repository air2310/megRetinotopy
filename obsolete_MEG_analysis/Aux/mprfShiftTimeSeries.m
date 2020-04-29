function y_sh = mprfShiftTimeSeries(y,sampling_rate, rad_or_samples, shift, F, do_plot)

if ~exist('do_plot','var') || isempty(do_plot)
    do_plot = false;
end


ft_y = fft(y,[],1);



if strcmpi(rad_or_samples,'radians')
    if ~exist('F','var') || isempty(F)
        error('Need to have a frequency')
    end
    
    shift = round((sampling_rate ./ F) .* (shift /( 2*pi)));
    shift = shift(:,ones(1,size(y,2)));
    
    
elseif strcmpi(rad_or_samples, 'samples')
    
    if shift ~=round(shift)
        warning('Shift is not a whole integer. Rounding')
        shift = round(shift);
    end
    
    if size(shift,2) == 1
        shift = shift(:,ones(1,size(y,2)));
        
    elseif shift(shift,2) == size(y,2)
        
        
        
    else
        error('Amount of shifts to do not match amount of time series')
        
    end
    
    
elseif strcmpi(rad_or_samples,'align_to_phase')
    if ~exist('F','var') || isempty(F)
        error('Need to have a frequency')
    end
    
    f_idx = mprfFreq2Index(size(y,1),F,sampling_rate);
    
    phase_diff = angle(ft_y(f_idx,:)) - shift;
    phase_diff(phase_diff > pi) = phase_diff(phase_diff > pi) - 2.*pi;
    phase_diff(phase_diff < -pi) = phase_diff(phase_diff < -pi) + 2.*pi;
    
    
    shift = round((sampling_rate ./ F) .* (phase_diff ./( 2.*pi)));
    
else
    error('Shift type not recognized');
    
    
end
    
Fs = 0:sampling_rate/size(y,1):sampling_rate/2; % Amount of repeates within
%Y for every frequency, given its length and sampling rate


ft_y = ft_y(1:floor(size(y,1)/2)+1,:);

ft_y_sh = ft_y.*exp(-1i.*2.*pi.*bsxfun(@times,Fs',shift).*(1./sampling_rate));

ft_y_sh = [ft_y_sh; conj(flipud(ft_y_sh(2:end-1,:)))];

y_sh = ifft(ft_y_sh,[],1,'symmetric');


if do_plot
    max_freq = 70;
    plot_idx = mprfFreq2Index(size(y,1),max_freq, sampling_rate);
    
    y_av = nanmean(y,2);
    y_av_sh = nanmean(y_sh,2);
    
    x_axis = linspace(0,max_freq,round(plot_idx));
    
    figure; hold on;
    plot(y_av,'k');
    plot(y_av_sh,'r');
    legend('Original','Shifted');
    
    f_idx = mprfFreq2Index(size(y,1),F,sampling_rate);
    
    
    y_ph = angle(ft_y(f_idx,:));
    y_ph_sh = angle(ft_y_sh(f_idx,:));
    
    figure; 
    plot(y_ph, y_ph_sh,'k.')
    axis([-pi pi -pi pi]);
    xlabel('Original phase')
    ylabel('Shifted phase')
    
    figure; hold on;
    plot(x_axis, log(nanmean(abs(ft_y(1:round(plot_idx),:)),2)),'k');
    plot(x_axis, log(nanmean(abs(ft_y_sh(1:round(plot_idx),:)),2)),'r');
    legend('Original','Shifted');
    xlabel('Frequency')
    ylabel('Log amplitude')
    
end


end
% 
% F = 34
% 
% sample_rate = 1000;
% phase_diff = pi;
% (2*pi) == (length(aaa) ./ sample_rate ./ F)
% 
% sample_diff = round((sample_rate ./ F) .* (phase_diff /( 2*pi)))
% 
% ft_a = fft(aaa);
% ph_a = angle(ft_a(round(mprfFreq2Index(length(aaa),F,sample_rate))));
% 
% bb = mprfShiftTimeSeries(aaa',sample_rate,sample_diff-1);
% 
% ft_b = fft(bb);
% ph_b = angle(ft_b(round(mprfFreq2Index(length(bb),F,sample_rate))));
% 
% ph_diff = ph_a - ph_b;
% ph_diff(ph_diff > pi) = ph_diff(ph_diff > pi) - 2*pi;
% ph_diff(ph_diff < -pi) = ph_diff(ph_diff < -pi) + 2*pi;
% 
% 
% disp(ph_diff)

% 
% 
% F = 1:551;
% 
% shift = (0.5 * 180) ./ (2.*pi .* F);
% 
% shift = 1;
% ft_aa = fft(aa);
% ph_a = angle(ft_aa(round(mprfFreq2Index(length(aa), F, 1000))));
% 
% bb = mprfShiftTimeSeries(aa',1000,shift);
% 
% ft_bb = fft(bb');
% ph_b = angle(ft_bb(round(mprfFreq2Index(length(bb), F, 1000))));
% 
% rad_diff_02 = ph_a - ph_b;
% 
% rad_diff = pi / (length(aa) / F);              % Radians per sample
% 
% disp(rad_diff)
% 
% th = (2*pi) .* F .* shift;
% 
% 
% return
% 
% 
% F = 11;
% 
% 
% ft_aa = fft(aa);
% ft_cc = fft(cc);
% 
% 
% ph_a = angle(ft_aa(round(mprfFreq2Index(length(aa), F, 1000))));
% ph_c = angle(ft_cc(round(mprfFreq2Index(length(cc), F, 1000))));
% 
% ph_diff_in = -(ph_a - ph_c);
% 
% rad_diff = pi / (length(aa)/2 / F)               % Radians per sample
% samp_diff = round(ph_diff_in .* (1/(pi / (length(aa)/2 / F))))      % Samples per radian
% 
% bb = mprfShiftTimeSeries(aa',1000, samp_diff/2);
% 
% ft_bb = fft(bb);
% ph_b = angle(ft_bb(round(mprfFreq2Index(length(bb), F, 1000))));
% 
% ft_bb = fft(bb);
% 
% ph_diff_out = (ph_a - ph_b);
% 
% ph_diff_out(ph_diff_out < -pi) = ph_diff_out(ph_diff_out < pi) + 2*pi;
% ph_diff_out(ph_diff_out > pi) = ph_diff_out(ph_diff_out > pi) - 2*pi;
% 
% disp(ph_diff_in);
% disp(ph_diff_out);
% 
% 
% 
% 
% 
% disp(ph_diff_in / (abs(ph_diff_out - ph_diff_in)));
% 
% 
% % 
% % t=0:0.001:0.1-0.001;        % Time
% % 
% % Fs = 1e3;                   % Sampling rate, i.e 1/step_size
% % 
% % freq1 = 100;
% %     
% % x1=cos(2*pi*freq1*t);       % Signal
% % 
% % Delay=10;                    % Delay
% % 
% % yp = fft(x1);               % FT of signal
% % 
% % yp = yp(1:length(x1)/2+1);  % Pruned
% % 
% % f = 0:Fs/length(x1):500;    % Frequencies of spectrum up to SR/2
% % 
% % yp = yp.*exp(-1i*2*pi*f*Delay*(1/Fs));
% % 
% % yp = [yp conj(fliplr(yp(2:end-1)))];
% % 
% % y = ifft(yp,'symmetric');
% % 
% % plot(t(1:100),x1(1:100),'b');
% % 
% % hold on;
% % 
% % plot(t(1:100),y(1:100),'r');
% % 
% % 
% % return
% % 
% % y = randn(1,10000);
% % sr_01 = 1000; % Assumed to be in Hz?
% % 
% % Freqs = 0:sr_01/length(y):sr_01/2; % Frequencies -> can't be right?? Only if sampling rate is in Hz
% % % 'Absolute' frequencies?? Not the repetitions per chunk of signal (y)....?
% % 
% % rad_D = -pi/4;
% % 
% % D = round(rad_D / (2*pi) * length(y));                   % Amount of points to shift
% % 
% % ft_y = fft(y);               % FT of signal
% % 
% % ft_y = ft_y(1:floor(length(y)/2)+1);  % Pruned
% % 
% % ft_y = ft_y.*exp(-1i.*2.*pi.*Freqs.*D.*(1./sr_01));
% % 
% % ft_y = [ft_y conj(fliplr(ft_y(2:end-1)))];
% % 
% % yy = ifft(ft_y,'symmetric');
% % 
% % figure;
% % plot(y,'b');
% % 
% % hold on;
% % 
% % plot(yy(1:length(y)),'r');
% % 
% 
% 
