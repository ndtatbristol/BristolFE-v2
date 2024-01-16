function mod = fn_tidy_L_model(mod)
%Removes unused nodes and elements but makes sure node numbers for relevant
%fields are correct
[mod.nds, mod.els, old_nds, ~] = fn_remove_unused_nodes(mod.nds, mod.els);
mod.bdry_lys2 = mod.bdry_lys2(old_nds);
mod.main_nd_i = mod.main_nd_i(old_nds);
mod.main_el_i = mod.main_el_i(old_nds);
end