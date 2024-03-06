function mod = fn_add_scatterer(mod, matls, scat_pts, scat_matl)
%SUMMARY
%   Adds scatterer to existing model by turning all elements inside
%   scat_pts to either matl(scat_matl) or void if = scat_matl


interface_el_name = 'ASI2D2';

%Remove interface elements if there are any
els_in_use = ~strcmp(mod.el_typ_i, interface_el_name);
[~, ~, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(els_in_use, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);


if scat_matl > 0
    mod = rmfield(mod, 'el_typ_i');
    mod = fn_set_els_inside_bdry_to_mat(mod, scat_pts, scat_matl);
else
    [~, els_in_use] = fn_elements_in_region2(mod.nds, mod.els, scat_pts);
    [~, ~, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(els_in_use, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);
    [mod.nds, mod.els, old_nds] = fn_remove_unused_nodes(mod.nds, mod.els);
    mod.bdry_lyrs = mod.bdry_lyrs(old_nds);
    mod.main_nd_i = mod.main_nd_i(old_nds);
end

%Add interface elements if needed
mod = fn_add_fluid_solid_interface_els(mod, matls);

%Set flag on which elements are within domain
mod.int_el_i = fn_elements_in_region(mod, mod.inner_bndry_pts);

end