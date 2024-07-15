function N = fn_symbolic_shape_function_matrix(n, no_dims, no_nds)
%n is no_nodes x 1 symbolic vector of shape functions; no_dims and no_nodes
%are scalar

N = sym('N', [no_dims no_dims * no_nds], 'real'); 
for i = 1:no_nds
    N(:, (i-1) * no_dims + 1: i * no_dims) = diag(repmat(n(i), [1, no_dims]));
end

end