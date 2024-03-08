function mod = fn_merge_subdomain_into_main_model(mainmod, submod, matls)

%Copy the mainmod mesh to the output
mod = mainmod;

%METHOD TOOD
%   Define bdry_nd_i = find(submod.bdry_lyrs == 1)
%   Create tmp_bdry_pts = submod.nds(bdry_nd_i, :)
%   Order tmp_bdry_pts so that they are in order to describe closed polygon
%   (should be function for this buried in the subdomain creation function)
%   Use this as the boundary below to delete elements from the two models
%   to guarantee consistency
%   Increase node reference in submod.els by n = size(mod.nds, 1)
%   Append remaining submod.nds to main model
%   Replace all references to original submod node numbers on boundary to
%   main_nd_i - this will take a but of thought to not mess up!



%main.doms{1}.mod.main_nd_i
%Remove elements outside inner boundary in submod
[not_in_use, ~] = fn_elements_in_region(submod, submod.inner_bndry_pts);
[old_els, new_els, submod.els, submod.el_mat_i, submod.el_abs_i, submod.el_typ_i] = fn_remove_unused_elements(not_in_use, submod.els, submod.el_mat_i, submod.el_abs_i, submod.el_typ_i);
% [submod.nds, submod.els, old_nds, new_nds] = fn_remove_unused_nodes(submod.nds, submod.els);

%Remove elements outside inner boundary from mod
[~, not_in_use] = fn_elements_in_region(mod, submod.inner_bndry_pts);
[old_els, new_els, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i] = fn_remove_unused_elements(not_in_use, mod.els, mod.el_mat_i, mod.el_abs_i, mod.el_typ_i);

keyboard
end
%combine the nodes
% nd_offset = size(submod.nds, 1);
% bdry_nds = find(mod.bdry_lys2 == 1) + nd_offset;
% mod.els = mod.els + nd_offset;
% submod.nds = [submod.nds; mod.nds];
% submod.els = [submod.els; mod.els];
% submod.el_mat_i = [submod.el_mat_i; mod.el_mat_i];
% submod.el_abs_i = [submod.el_abs_i; mod.el_abs_i];
% submod.el_typ_i = [submod.el_typ_i; mod.el_typ_i];
% 
% %they are now combined, but not actually joined. Identify duplicate
% %nodes and remove
% duplicates = [];
% for i = 1:numel(bdry_nds)
%     %find all nodes within a reasonable distance, then check these for
%     %duplicates
%     tmp = find(sqrt(sum((submod.nds - submod.nds(bdry_nds(i), :)) .^ 2, 2)) < submod.el_size * 2);
%     % plot(val_mod.nds(tmp, 1), val_mod.nds(tmp, 2), 'g.');
%     for j = 1:numel(tmp)
%         tmp2 = find(sqrt(sum((submod.nds - submod.nds(tmp(j), :)) .^ 2, 2)) < submod.el_size / 10)';
%         if numel(tmp2) == 2
%             % plot(val_mod.nds(tmp(j), 1), val_mod.nds(tmp(j), 2), 'co');
%             duplicates = [duplicates; tmp2];
%         end
%     end
% end
% for i = 1:size(duplicates, 1)
%     submod.els(submod.els == duplicates(i, 1)) = duplicates(i, 2);
% end