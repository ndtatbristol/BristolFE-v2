function [history_output, field_output, force_output, field_output_time] = fn_explicit_dynamic_solver_v5(...
    K, C, M, time, ...
    forcing_indices, forcing_functions, ...
    disp_indices, disp_functions, ...
    history_indices, field_output_every_n_frames, varargin)
%SUMMARY
%   Solves explicit dynamic FE problem given applied displacements or
%   applied forces
%INPUTS
%   K - m x m global stiffness matrix
%   C - m x m global damping matrix
%   M - m x m global mass matrix
%   time = n-element vector of times
%   forcing_indices - p-element vector of global matrix indices at which 
%   forcing_functions will be applied (i.e. force input)
%   forcing_functions - p x n matrix of forces to apply
%   disp_indices - r-element vector of global matrix indices at which 
%   disp_functions will be applied (i.e. displacement input)
%   disp_functions - r x n matrix of displacments to apply
%   history_indices - q-element vector of global matrix indices for complete 
%   time-history outputs
%   field_output_every_n_frames - complete displacement field will be
%   output every n frames (set to inf for no field output
%OUTPUTS
%   history_output - q x n matrix of time histories
%   field_output - m x floor(n / field_output_every_n_frames) matrix of displacements at all nodes
%   force_output - r x n matrix for force histories at points where
%   displacements are imposed (empty if no displacement input is used)

%--------------------------------------------------------------------------
field_output_is_vel = 1;
if isempty(varargin)
	use_gpu_if_present = 1;
else
	use_gpu_if_present = varargin{1};
end

gpu_present = fn_test_if_gpu_present_and_working;
if use_gpu_if_present && gpu_present
	use_gpu = 1;
    reset(gpuDevice);
else
	use_gpu = 0;
end

%Error checks
if size(K,1) ~= size(K, 2)
    error('K must be square matrix');
end
if size(C,1) ~= size(C, 2)
    error('C must be square matrix');
end
if size(M,1) ~= size(M, 2)
    error('M must be square matrix');
end

fprintf('Explicit time marching v5 (GPU = %i, time steps = %i) ', use_gpu,  numel(time));
time_step = time(2) - time(1);

%initialise history and field output variables
if isempty(history_indices)
    history_output = [];
else
    history_output = zeros(length(history_indices), length(time));
end

if isempty(disp_indices)
    force_output = [];
else
    force_output = zeros(length(disp_indices), length(time));
    tmp = disp_functions;
    tmp = [zeros(size(disp_functions, 1), 2), disp_functions];
    accn = zeros(size(disp_functions));
    accn = (tmp(:, 3:end) - 2 * tmp(:, 2:end-1) + tmp(:, 1:end-2)) / time_step ^ 2;
end

if ~isinf(field_output_every_n_frames)
    field_output_ti = 1:field_output_every_n_frames:length(time);
    field_output = zeros(size(K, 1), length(field_output_ti));
    field_output_time = zeros(1, length(field_output_ti));
else
    field_output_ti = [];
    field_output = [];
    field_output_time = [];
end

diag_M = spdiags(sum(M).', 0, size(K,1), size(K,2));
inv_M = spdiags(1 ./ sum(M).', 0, size(K,1), size(K,2));

u_previous = zeros(size(K, 1), 1);
u_dot_previous = zeros(size(K, 1), 1);
u_dot_dot_previous = zeros(size(K, 1), 1);

f = zeros(size(K, 1), 1);

if use_gpu
	K = gpuArray(K);
	C = gpuArray(C);
    inv_M = gpuArray(inv_M);
    diag_M = gpuArray(diag_M);
	u_previous = gpuArray(u_previous);
	u_dot_previous = gpuArray(u_dot_previous);
    u_dot_dot_previous = gpuArray(u_dot_dot_previous);
    f = gpuArray(f);
    if ~isempty(forcing_indices)
        forcing_functions = gpuArray(forcing_functions);
        forcing_indices = gpuArray(forcing_indices);
    end
    if ~isempty(disp_indices)
        disp_indices = gpuArray(disp_indices);
        force_output = gpuArray(force_output);
    end
	if ~isempty(history_indices)
		history_indices = gpuArray(history_indices);
		history_output = gpuArray(history_output);
    end
	if ~isinf(field_output_every_n_frames)
		field_output = gpuArray(field_output);
	end
end

%Main time marching loop
t1 = clock;
ti_start = inf;
if ~isempty(forcing_indices)
    ti_start = min(min(find(sum(forcing_functions))), ti_start);
end
if ~isempty(disp_indices)
    ti_start = min(min(find(sum(disp_functions))), ti_start);
end
for ti = ti_start:length(time)
    %set force at forcing node equal to excitation signal at this instant
    %in time
    if ti <= size(forcing_functions, 2)
        f(forcing_indices) = forcing_functions(:, ti);
    end

    %work out acceleration
    u_dot_dot_previous = inv_M * (f - K * u_previous - C * u_dot_previous);

    %work out velocity
    u_dot = u_dot_previous + time_step * u_dot_dot_previous;
    
    %work out displacement at next time step
    u = u_previous + time_step * u_dot;
    
    %impose displacements
    if ~isempty(disp_indices)
        u(disp_indices) = disp_functions(:, ti);
        force_output(:, ti) = diag_M(disp_indices, disp_indices) * accn(:, ti) + K(disp_indices, :) * u_previous;
    end

    %history output
    if ~isempty(history_indices)
        history_output(:, ti) = u(history_indices);
    end
    
    %field output
    [tmp, fi] = ismember(ti, field_output_ti);
    if tmp
        field_output_time(fi) = time(ti);
        if field_output_is_vel
            field_output(:, fi) = u_dot;
        else
            field_output(:, fi) = u;
        end
    end
    
    %overwrite previous values with current ones ready for next loop
    u_previous = u;
    u_dot_previous = u_dot;
    
    %Show how far through calculation is
    if rem(ti,100) == 0
        fn_simple_text_progress_bar(ti, length(time));
    end
end
if use_gpu
	if ~isempty(history_indices)
		history_output = gather(history_output);
	end
	if ~isempty(disp_indices)
		force_output = gather(force_output);
	end
	if ~isinf(field_output_every_n_frames)
		field_output = gather(field_output);
    end
    reset(gpuDevice);
end

t2 = etime(clock, t1);
fprintf('completed in %.2f secs\n', t2);

end