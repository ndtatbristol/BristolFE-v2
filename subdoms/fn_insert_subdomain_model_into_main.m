function mod = fn_insert_subdomain_model_into_main(mn_mod, dm_mod, matls)
%Can delete matls argument after debugging - used for plotting only
% figure;h = fn_show_geometry(dm_mod, matls, []);

mod = mn_mod;

%Remove the elements from model that are inside region
els_in_use = ones(size(mod.els, 1), 1);
els_in_use(dm_mod.main_int_el_i) = 0;
[~, ~, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(els_in_use, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);

% figure;h = fn_show_geometry(mod, matls, []);

%stick the subdomain nodes at the end of the main model nodes
dm_nd_offset = size(mod.nds, 1)
mod.nds = [mod.nds; dm_mod.nds];

%remove elements in subdomain that are outside region
[~, ~, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i] = fn_remove_unused_elements(dm_mod.int_el_i, dm_mod.els, dm_mod.el_mat_i, dm_mod.el_abs_i, dm_mod.el_typ_i);

% figure;h = fn_show_geometry(dm_mod, matls, []);

%Looks OK to here, but then coordinates of some/all nodes in subdomain go wrong

mod.els =      [mod.els;      dm_mod.els + dm_nd_offset];
mod.el_mat_i = [mod.el_mat_i; dm_mod.el_mat_i];
mod.el_abs_i = [mod.el_abs_i; dm_mod.el_abs_i];
mod.el_typ_i = [mod.el_typ_i; dm_mod.el_typ_i];

figure;h = fn_show_geometry(mod, matls, []);
end

