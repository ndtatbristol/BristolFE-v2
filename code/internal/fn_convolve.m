function h = fn_convolve(f, g, dim, varargin)

if numel(varargin) < 1
    use_gpu = 1;
else
    use_gpu = varargin{1};
end

%Convolves f by g along specified dim
fprintf('    Convolving %i signals with %i pts (GPU = %i)', prod(size(f)) / size(f,dim), size(f,dim), use_gpu);
t1 = clock;

fft_pts = size(f, dim) * 2;

if use_gpu
    f = gpuArray(f);
    g = gpuArray(g);
end

F = fft(f, fft_pts, dim);
G = fft(g, fft_pts, dim);
H = F .* G;

h = ifft(H, fft_pts, dim);
h = h(1: size(f,1), 1: size(f,2));

if use_gpu
    h = gather(h);
    reset(gpuDevice);
end

t2 = etime(clock, t1);
fprintf(' .......... completed in %.2f secs\n', t2);
end