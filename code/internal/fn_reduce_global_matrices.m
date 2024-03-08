function [global_matrix_nodes, global_matrix_dofs, Q, varargout] = fn_reduce_global_matrices(global_matrix_nodes, global_matrix_dofs, Q, varargin)

%Reduces global matrices by knocking out all empty rows and cols. The
%lookup vectors, global_matrix_nodes and global_matrix_dofs, are also
%updated. Mats is a cell array, typically containing M, K and possibly C.

for m = 1:length(varargin)
    if m == 1
        valid_indices = sum(abs(varargin{m})) > 0;
    else
        valid_indices = valid_indices | sum(abs(varargin{m})) > 0;
    end
end
for m = 1:length(varargin)
    varargout{m} = varargin{m}(valid_indices, valid_indices);
end
if ~isempty(Q) %needed for legacy code that used Q matrix
    Q = Q(:, valid_indices);
end
global_matrix_nodes = global_matrix_nodes(valid_indices);
global_matrix_dofs  = global_matrix_dofs(valid_indices);

end