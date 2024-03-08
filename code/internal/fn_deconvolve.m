function h = fn_deconvolve(f, g, dim)

%De-convolves f by g along specified dim

fft_pts = size(f, dim) * 2;
% fft_pts = size(f, dim);

F = fft(f, fft_pts, dim);
G = fft(g, fft_pts, dim);
H = F ./ G;

h = ifft(H, fft_pts, dim);
h = h(1: size(f,1), 1: size(f,2));
end