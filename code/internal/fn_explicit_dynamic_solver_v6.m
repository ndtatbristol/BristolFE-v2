function [history_output, field_output, force_output, field_output_time] = fn_explicit_dynamic_solver_v6(...
    K, C, M, time, ...
    forcing_indices, forcing_functions, ...
    disp_indices, disp_functions, ...
    history_indices, field_output_every_n_frames, varargin)
%v6 changes - splitting time stepping into two methods depending on whether
%the velocities are calculated on the preceding half step (old method) or
%on current step (should be more stable, but slower).

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
if numel(varargin) < 1
	use_gpu_if_present = 1;
else
	use_gpu_if_present = varargin{1};
end
if numel(varargin) < 2
	field_output_type = 'KE';
else
	field_output_type = varargin{2};
end
if numel(varargin) < 3
	solver_mode = 'vel at last half time step';
else
	solver_mode = varargin{3};
end

%--------------------------------------------------------------------------
switch field_output_type
    case {'mean(u1)', 'mean(u2)', 'mean(u3)', 'curl(u)', 'div(u)', 'raw(u)'}
        field_output_is_vel = 0;
    otherwise
        field_output_is_vel = 1;
end

gpu_present = fn_test_if_gpu_present_and_working;
if use_gpu_if_present && gpu_present
	use_gpu = 1;
    reset(gpuDevice);
else
	use_gpu = 0;
end

%Error checks
if size(K, 1) ~= size(K, 2)
    error('K must be square matrix');
end
if size(C,1) ~= size(C, 2)
    error('C must be square matrix');
end
if size(M,1) ~= size(M, 2)
    error('M must be square matrix');
end

ndf = size(K, 1);

fprintf('Explicit time marching v6 (GPU = %i, time steps = %i, DOF = %i) ', use_gpu,  numel(time), ndf);
dt = time(2) - time(1);

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
    accn = (tmp(:, 3:end) - 2 * tmp(:, 2:end-1) + tmp(:, 1:end-2)) / dt ^ 2;
end

if ~isinf(field_output_every_n_frames)
    field_output_ti = 1:field_output_every_n_frames:length(time);
    field_output = zeros(ndf, length(field_output_ti));
    field_output_time = zeros(1, length(field_output_ti));
else
    field_output_ti = [];
    field_output = [];
    field_output_time = [];
end

diag_M = spdiags(sum(M).', 0, ndf, ndf);
inv_M = spdiags(1 ./ sum(M).', 0, ndf, ndf);

u_minus_1 = zeros(ndf, 1);
u_minus_2 = zeros(ndf, 1);

f = zeros(ndf, 1);

switch solver_mode 
    case 'vel at last half time step'
        B1 = 2 * speye(ndf) - dt * inv_M * C - dt ^ 2 * inv_M * K;
        B2 = dt * inv_M * C - speye(ndf);
        B3 = 0;
    case 'vel at curent time step'
        B1 = 2 * speye(ndf) - dt ^ 2 * inv_M * K;
        B2 =    -speye(ndf) + dt / 2 * inv_M * C;
        B3 =     speye(ndf) + dt / 2 * inv_M * C;
end

if use_gpu
	% K = gpuArray(K);
	% C = gpuArray(C);
    inv_M = gpuArray(inv_M);
    diag_M = gpuArray(diag_M);
	u_minus_1 = gpuArray(u_minus_1);
	u_minus_2 = gpuArray(u_minus_2);
    B1 = gpuArray(B1);
    B2 = gpuArray(B2);
    f = gpuArray(f);
end

%Main time marching loop
t1 = clock;
ti_start = inf;
if ~isempty(forcing_indices)
    ti_start = min(min(find(sum(abs(forcing_functions)))), ti_start);
end
if ~isempty(disp_indices)
    ti_start = min(min(find(sum(abs(disp_functions)))), ti_start);
end

prog_dot_ti = interp1(linspace(0, 1, length(time) - ti_start + 1), ti_start:length(time), linspace(0,1,11), 'nearest');
prog_dot_ti = prog_dot_ti(2: end);


for ti = ti_start:length(time)
    %set force at forcing node equal to excitation signal at this instant
    %in time
    if ti <= size(forcing_functions, 2)
        f(forcing_indices) = forcing_functions(:, ti);
    end

    %Main calculation!
    switch solver_mode 
        case 'vel at last half time step'
            u =       dt ^ 2 * inv_M * f + B1 * u_minus_1 + B2 * u_minus_2;
        case 'vel at curent time step'
            u = B3 \ (dt ^ 2 * inv_M * f + B1 * u_minus_1 + B2 * u_minus_2);
    end
       
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
            field_output(:, fi) = (u - u_minus_1) / dt;
        else
            field_output(:, fi) = u;
        end
    end
    
    %overwrite previous values with current ones ready for next loop
    u_minus_2 = u_minus_1;
    u_minus_1 = u;
    
    %Show how far through calculation is
    if ismember(ti, prog_dot_ti)
        fprintf('.')
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
fprintf(' completed in %.2f secs\n', t2);

end