dt = 0.001;
t = dt:dt:1;
fs = (0:length(t)-1)/max(t);

bl_s = randn(size(t));
bb_s = randn(size(t));
sl_s = sin(2*pi*t*10);
bl_w = 10;
bb_w = 5;
sl_w = 15;

s1 = bl_w*bl_s + bb_w*bb_s + sl_w*sl_s;
s2 = cumsum(s1);

figure,
subplot(2,1,1)
plot(t,s1, t,s2);

subplot(2,1,2)
plot(fs,abs(fft(s1)), fs,abs(fft(s2)));
set(gca, 'XScale', 'log', 'YScale', 'log', 'XLim', [2 200])

