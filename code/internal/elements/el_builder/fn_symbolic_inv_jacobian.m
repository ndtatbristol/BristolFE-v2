function [detJ, invJ_times_detJ] = fn_symbolic_inv_jacobian(n, nds, q)
no_nds = size(nds, 1);
no_dims = size(nds, 2);
phys_dims = 3;

%Shape function matrix N for coordinates (may not be same size as for other
%quantities, hence why it is done separately here)
N = fn_symbolic_shape_function_matrix(n, no_dims, no_nds);

%Phys coordinates in terms of natural ones
x = N * reshape(nds', [], 1);

%Returns 3x3 inverse Jacobian regardless of dimensionality of inputs
J = jacobian(x, q);
detJ = det(J);
invJ = inv(J);
invJ_times_detJ = invJ * detJ;
% invJ_times_detJ = [ ...
%     invJ_times_detJ,                            sym(zeros(no_dims,phys_dims - no_dims)); ... 
%     sym(zeros(phys_dims - no_dims, no_dims)),   sym(eye(phys_dims - no_dims))];

end
