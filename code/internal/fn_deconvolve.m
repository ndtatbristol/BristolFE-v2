function h = fn_deconvolve(f, g, dim, varargin)
%De-convolves f by g along specified dim

if numel(varargin) < 1
    use_gpu_if_present = 1;
else
    use_gpu_if_present = varargin{1};
end

use_gpu = 0;
if use_gpu_if_present
    gpu_present = fn_test_if_gpu_present_and_working;
    if gpu_present
	    use_gpu = 1;
        reset(gpuDevice);
    end
end

fprintf('    Deconvolving %i signals with %i pts (GPU = %i)', prod(size(f)) / size(f,dim), size(f,dim), use_gpu);
t1 = clock;

fft_pts = size(f, dim) * 2;

if use_gpu
    f = gpuArray(f);
    g = gpuArray(g);
end


F = fft(f, fft_pts, dim);
G = fft(g, fft_pts, dim);
H = F ./ G;

h = ifft(H, fft_pts, dim);
h = h(1: size(f,1), 1: size(f,2));

if use_gpu
    h = gather(h);
    reset(gpuDevice);
end

t2 = etime(clock, t1);
fprintf(' .......... completed in %.2f secs\n', t2);

end