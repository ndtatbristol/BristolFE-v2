function gl_ind = fn_gl_ind_for_nd_and_dof(gl_lookup, nds, dofs)
gl_ind = gl_lookup(sub2ind(size(gl_lookup),nds, dofs));
end