function mod = fn_insert_subdomain_model_into_main(mn_mod, dm_mod, matls)
%Can delete matls argument after debugging - used for plotting only
% figure;h = fn_show_geometry(dm_mod, matls, []);

if isfield(mn_mod, 'el_abs_i')
    mn_mod.el_abs_i = mn_mod.el_abs_i;
else
    mn_mod.el_abs_i = zeros(size(mn_mod.el_mat_i));
end

mod = mn_mod;



%Remove the elements from mmain odel that are inside region
els_in_use = ones(size(mod.els, 1), 1);
els_in_use(dm_mod.main_int_el_i) = 0;
[~, ~, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(els_in_use, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);

% figure;h = fn_show_geometry(mod, matls, []);

%Stick the subdomain nodes at the end of the main model nodes and remember
%the offset that needs to be added to references to these
dm_nd_offset = size(mod.nds, 1);
mod.nds = [mod.nds; dm_mod.nds];


%remove elements in subdomain that are outside region
[~, ~, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i] = fn_remove_unused_elements(dm_mod.int_el_i, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i);

% figure;h = fn_show_geometry(dm_mod, matls, []);

%Change boundary node references in sub-domain model to main node
%equivalents
ed_nds = find(dm_mod.bdry_lyrs == 1);
i = ismember(dm_mod.els, ed_nds);
dm_mod.els(i) = interp1(1:size(dm_mod.nds, 1), dm_mod.main_nd_i, dm_mod.els(i), 'nearest') - dm_nd_offset;
dm_mod.els(dm_mod.els == 0) = -dm_nd_offset; %need to handle zeros (which need to be still zeros when offset)
mod.els =      [mod.els;      dm_mod.els + dm_nd_offset];
mod.el_mat_i = [mod.el_mat_i; dm_mod.el_mat_i];
mod.el_abs_i = [mod.el_abs_i; dm_mod.el_abs_i];
mod.el_typ_i = [mod.el_typ_i; dm_mod.el_typ_i];

% figure;h = fn_show_geometry(mod, matls, []);
end


