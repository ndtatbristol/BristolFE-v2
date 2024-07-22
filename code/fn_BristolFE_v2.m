function [res, mats] = fn_BristolFE_v2(mod, matls, steps, fe_options)
%SUMMARY
%   Entry function for Bristol FE v2.
%INPUTS
%   mod - description of mesh including nodes, elements, material
%   indices, and possibly absorbing indices if absorbing layers are used.
%   matls - description of materials used in mod
%   steps - description of one or more (use cell array) steps in which 
%       loads are applied, including details of the load and what is
%       recorded
%OUTPUTS
%   res - results from each load step
%   mats - global matrices for model
%--------------------------------------------------------------------------

default_options.use_gpu_if_available = 1;
default_options.field_output_every_n_frames = inf;
default_options.global_matrix_builder_version = 'v5';
default_options.dynamic_solver_version = 'v6';
default_options.solver_mode = 'vel at last half time step';
default_options.field_output_type = 'KE';
fe_options = fn_set_default_fields(fe_options, default_options);

%Check inputs - are all mesh and material details consistent?

%If mod.el_typ_i not specified at this point, generate it based on matls
if ~isfield(mod, 'el_typ_i')
    mod.el_typ_i = {matls(mod.el_mat_i).el_typ};
    mod.el_typ_i = mod.el_typ_i(:);
end

if ~isfield(mod, 'el_abs_i')
    mod.el_abs_i = zeros(size(mod.el_mat_i));
end


%Build the global matrices
switch fe_options.global_matrix_builder_version
    case 'v4'
        [mats.K, mats.C, mats.M, mats.gl_lookup] = fn_build_global_matrices_v4(mod.nds, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i, matls, fe_options);
    otherwise %v5 is default
        [mats.K, mats.C, mats.M, mats.gl_lookup] = fn_build_global_matrices_v5(mod.nds, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i, matls, fe_options);
end

if isempty(steps)
    %Useful if only matrices are required
    res = {};
    return
end

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
            [frc_gi, ~, ~, valid_appl_frcs] = fn_nds_and_dfs_to_gi(steps{s}.load.frc_nds, steps{s}.load.frc_dfs, mats.gl_lookup);
            frcs = steps{s}.load.frcs;
            if isfield(steps{1}.load, 'wts')
                frcs = steps{1}.load.wts(:) * frcs;
            end
            if size(frcs, 1) > 1
                frcs = frcs(valid_appl_frcs, :);
            end
        else
            frc_gi = [];
            frcs = [];
        end
        if isfield(steps{s}.load, 'dsps')
            %[dsp_gi, res{s}.frc_nds, res{s}.frc_dfs, res{s}.valid_appl_dsps] = fn_nds_and_dfs_to_gi(steps{s}.load.dsp_nds, steps{s}.load.dsp_dfs, mats.gl_lookup);
            [dsp_gi, ~, ~, res{s}.valid_appl_dsps] = fn_nds_and_dfs_to_gi(steps{s}.load.dsp_nds, steps{s}.load.dsp_dfs, mats.gl_lookup);
            dsps = steps{s}.load.dsps;
            if isfield(steps{1}.load, 'wts')
                dsps = steps{1}.load.wts(:) * dsps;
            end
            if size(dsps,1) > 1
                dsps = dsps(res{s}.valid_appl_dsps, :);
            end
        else
            dsp_gi = [];
            dsps = [];
        end
        if isfield(steps{s}.mon, 'nds')
            % [hist_gi, res{s}.dsp_nds, res{s}.dsp_dfs, res{s}.valid_mon_dsps] = fn_nds_and_dfs_to_gi(steps{s}.mon.nds, steps{s}.mon.dfs, mats.gl_lookup);
            [hist_gi, ~, ~, res{s}.valid_mon_dsps] = fn_nds_and_dfs_to_gi(steps{s}.mon.nds, steps{s}.mon.dfs, mats.gl_lookup);
        else
            hist_gi = [];
        end

        %Explicit dynamic analysis - (2) run it!
        if numel(steps) > 1
            fprintf('    (%2i/%2i) ', s, numel(steps));
        else
            fprintf('    ');
        end
        switch fe_options.dynamic_solver_version
            case'v5'
                [mon_dsps, fld, mon_frcs, res{s}.fld_time] = ...
                    fn_explicit_dynamic_solver_v5(mats.K, mats.C, mats.M, t, ...
                    frc_gi, frcs, dsp_gi, dsps, hist_gi, fe_options.field_output_every_n_frames, fe_options.use_gpu_if_available);
            otherwise %v6 is now the default
                [mon_dsps, fld, mon_frcs, res{s}.fld_time] = ...
                    fn_explicit_dynamic_solver_v6(mats.K, mats.C, mats.M, t, ...
                    frc_gi, frcs, dsp_gi, dsps, hist_gi, fe_options.field_output_every_n_frames, fe_options.use_gpu_if_available, fe_options.field_output_type, fe_options.solver_mode);
        end
        
        %Parse the monitored history outputs
        if isfield(steps{s}.mon, 'nds')
            res{s}.dsps(res{s}.valid_mon_dsps, :) = mon_dsps;
            res{s}.dsp_gi = hist_gi; %note that this is NOT dsp_gi above as that is indices of applied displacements (confusing!)
        else
            res{s}.dsps = [];
            res{s}.dsp_gi = [];
        end
        if isfield(steps{s}.load, 'dsps')
            res{s}.frcs(res{s}.valid_appl_dsps, :) = mon_frcs;
            res{s}.frc_gi = frc_gi;
        else
            res{s}.frcs = [];
            res{s}.frc_gi = [];
        end

        %Convert field output to element values
        if ~isempty(fld)
            % res{s}.fld = fn_get_plot_vals_v3(fld, mats.gl_lookup, mod.els, mats.M);
            % res{s}.fld = fn_get_plot_vals_v4(fld, mats.gl_lookup, mod.nds, mod.els, mats.M, fe_options.field_output);
            res{s}.fld = fn_get_field_output(fld, mod, mats, fe_options.field_output_type);
        else
            res{s}.fld = [];
        end
        
    end
end

end

function [gi, nds, dfs, valid] = fn_nds_and_dfs_to_gi(nds, dfs, gl_lookup)
gi = zeros(numel(nds), 1);
for i = 1:numel(nds)
    gi(i) = gl_lookup(nds(i), dfs(i));
end
valid = gi > 0;
gi = gi(valid);
nds = nds(valid);
dfs = dfs(valid);
end


