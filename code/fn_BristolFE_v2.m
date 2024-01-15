function [res, mats] = fn_BristolFE_v2(mod, matls, steps, options)
%SUMMARY
%   Entry function for Bristol FE v2.
%INPUTS
%   mod - description of mesh including nodes and elements
%   matls - description of materials used in mod
%   steps - description of one or more (use cell array) steps in which 
%       loads are applied, including details of the load and what is
%       recorded
%OUTPUTS
%   res - results from each load step
%   mats - global matrices for model
%--------------------------------------------------------------------------

default_options.use_gpu_if_present = 1;

options = fn_set_default_fields(options, default_options);

%Check inputs - are all mesh and material details consistent?


%Build the global matrices
[mats.K, mats.C, mats.M, mats.gl_lookup] = fn_build_global_matrices_v4(mod.nds, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i, matls, options);

%Convert steps into cell array if not already
if ~iscell(steps)
    steps = {steps};
end

%Loop over each step, determine type of analysis to do and then do it!
for s = 1:numel(steps)
    if isfield(steps{s}.load, 'time')
        %Explicit dynamic analysis - (1) sort out inputs
        t = steps{s}.load.time;
        if isfield(steps{s}.load, 'frcs')
            frc_gi = fn_nds_and_dfs_to_gi(steps{s}.load.frc_nds, steps{s}.load.frc_dfs, mats.gl_lookup);
            frcs = steps{s}.load.frcs;
        else
            frc_gi = [];
            frcs = [];
        end
        if isfield(steps{s}.load, 'dsps')
            dsp_gi = fn_nds_and_dfs_to_gi(steps{s}.load.dsp_nds, steps{s}.load.dsp_dfs, mats.gl_lookup);
            dsps = steps{s}.load.dsps;
        else
            dsp_gi = [];
            dsps = [];
        end
        if isfield(steps{s}.mon, 'nds')
            hist_gi = fn_nds_and_dfs_to_gi(steps{s}.mon.nds, steps{s}.mon.dfs, mats.gl_lookup);
        else
            hist_gi = [];
        end
        if isfield(steps{s}.mon, 'field_output_every_n_frames')
            f_every_n = steps{s}.mon.field_output_every_n_frames;
        else
            f_every_n = inf;
        end

        %Explicit dynamic analysis - (2) run it!
        [res{s}.dsp, fld, res{s}.frc, res{s}.fld_time] = fn_explicit_dynamic_solver_v5(mats.K, mats.C, mats.M, t, frc_gi, frcs, dsp_gi, dsps, hist_gi, f_every_n, options.use_gpu_if_present);
        
        %Convert field output to element values
        res{s}.fld = fn_get_plot_vals_v3(fld, mats.gl_lookup, mod.els, mats.M);
        
    end
end

end

function gi = fn_nds_and_dfs_to_gi(nds, dfs, gl_lookup)
gi = zeros(numel(nds), 1);
for i = 1:numel(nds)
    gi(i) = gl_lookup(nds(i), dfs(i));
end
gi = gi(gi > 0);
end


