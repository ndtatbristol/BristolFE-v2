function [hist_gi, doms] = fn_find_gl_inds_for_domains(gl_lookup, doms)
hist_gi = [];
for d = 1:numel(doms)
    for s = 1:numel(doms{d}.scats)
        k = find(doms{d}.scats{s}.mod.bdry_lys2 > 0); %these are the domain nodes on the boundary layers where data is needed
        [tmp, ~, ~, nd_i] = fn_global_indices_for_all_dof_at_nodes(gl_lookup, doms{d}.scats{s}.mod.main_nd_i(k)); %these are associate global indices for all DOF at these nodes
        doms{d}.scats{s}.mod.main_hist_gi = [numel(hist_gi) + 1: numel(hist_gi) + numel(tmp)]'; %this is where to look up boundary data for this domain in the history output of main model
        doms{d}.scats{s}.mod.bdry_lys_gi = doms{d}.scats{s}.mod.bdry_lys2(k(nd_i)); %these are the associated layer indices
        hist_gi = [hist_gi; tmp]; %this accumulates main model global DOF where data is required
    end
end
end