function gl_lookup = fn_create_fast_lookup(global_matrix_nodes, global_matrix_dofs, max_nds, max_dofs)
%returns a matrix of global matrix lookup indices where rows are nodes and
%cols are dofs
ind = 1:numel(global_matrix_nodes);
gl_lookup = zeros(max(max(global_matrix_nodes(:)), max_nds), max(max(global_matrix_dofs(:)), max_dofs));
for i = 1:length(ind)
    gl_lookup(global_matrix_nodes(i), global_matrix_dofs(i)) = ind(i);
end
end