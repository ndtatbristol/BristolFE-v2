function h = fn_convolve(f, g, dim, varargin)
%Convolves f by g along specified dim

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

fn_console_output(sprintf('Convolving %i signals with %i pts (GPU = %i) ... ', prod(size(f)) / size(f,dim), size(f,dim), use_gpu));
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

fn_console_output(sprintf('completed in %.2f secs\n', etime(clock, t1)), [], 0);
end