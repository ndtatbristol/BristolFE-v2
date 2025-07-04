function varargout = fn_FE_entry_point(mod, matls, steps, fe_options)
%SUMMARY
%   Common entry point for different FE solvers.
%USAGE
%   res = fn_FE_entry_point(mod, matls, steps, fe_options)
%   [res, mats] = fn_FE_entry_point(mod, matls, steps, fe_options)
%   fe_options = fn_FE_entry_point([], [], [], fe_options)
%INPUTS
%   mod - description of mesh including nodes, elements, material
%   indices, and possibly absorbing indices if absorbing layers are used.
%   matls - description of materials used in mod
%   steps - description of one or more (use cell array) steps in which
%       loads are applied, including details of the load and what is
%       recorded
%OUTPUTS
%   res - results from each load step
%   [mats - global matrices for model]
%   fe_options - special case used to obtain options (including defaults)

%--------------------------------------------------------------------------
%FE_OPTIONS meanings and defaults
%Default solver
default_options.solver = 'BristolFE';
%Default properties for absorbing layers (used to map fractional
%distance into absorbing layer in range 0 to 1 into modifications to
%elements stiffness and damping matrices
default_options.damping_power_law = 3;
default_options.max_damping = 3.1415e+07;
default_options.max_stiffness_reduction = 0.01;
%Solver precision
default_options.solver_precision = 'double';
%How often to output field output (inf = never)
default_options.field_output_every_n_frames = inf;
%Which DoF to include in model, use [] for all
default_options.dof_to_use = []; 

%--------------------------------------------------------------------------
fe_options = fn_set_default_fields(fe_options, default_options);
global COMMENT_INDENT_LEVEL 

%Set the solver
switch fe_options.solver
    case 'Abaqus'
        fn_solver = @fn_ABAQUS;
    case 'BristolFE'
        fn_solver = @fn_BristolFE_v2;
    case 'pogo'
        fn_solver = @fn_pogoFE;
    otherwise
        error('Invalid FE solver');
end

%Call selected solver with appropriate number of outputs
t1 = clock;
fn_console_output(['Starting FE solver (', fe_options.solver, ')\n'])
fn_increment_indent_level;
if nargout == 1
    varargout{1} = fn_solver(mod, matls, steps, fe_options);
elseif nargout == 2
    [varargout{1}, varargout{2}] = fn_solver(mod, matls, steps, fe_options);
end
fn_decrement_indent_level;
fn_console_output(sprintf(['FE solver (', fe_options.solver, ') completed in %.2f secs\n'], etime(clock, t1)));
end



