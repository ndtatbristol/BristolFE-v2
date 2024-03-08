function h = fn_convolve(f, g, dim)

%Convolves f by g along specified dim

fft_pts = size(f, dim) * 2;

F = fft(f, fft_pts, dim);
G = fft(g, fft_pts, dim);
H = F .* G;

h = ifft(H, fft_pts, dim);
h = h(1: size(f,1), 1: size(f,2));
end