function s = fn_gaussian_pulse(t, centre_freq, no_cycles)
%SUMMARY
%   Returns a Gaussian pulse with specifed centre frequency and number of
%   cycles (based on -40dB points) given a specified time axis.

sz= size(t);
t = t(:);
beta = 0.01;
ct = mean(t);
T = no_cycles / (2 * centre_freq * sqrt(log(1/beta)));
env = exp(-((t - ct)/T) .^ 2);
i = min(find(env > 0.01));
s = env .* sin(2*pi * centre_freq * (t - ct));
s = [s(i:end); zeros(i - 1, 1)];
s = reshape(s, sz);
end