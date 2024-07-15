function [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_ASI2D2(nds, els, D, rho, varargin)
%SUMMARY
%	This function was created automatically by fn_create_element_matrix_file
%	and contains code to return the stiffness and mass matrices
%	for multiple elements of the same material and type given by the latter
%	part of the filename, fn_el_ASI2D2.
%INPUTS
%	nds - n_nds x n_dims matrix of nodal coordinates
%	els - n_els x n_nds_per_el matrix of node indices for each elements
%	D - ns x ns material stiffness matrix
%	rho - material density
%	[dofs_to_use = [] - optional string listing the DoFs to use, e.g. '12'. Use [] for all]
%OUTPUTS
%	el_K, el_C, el_M - n_els x n_dfs_per_el x n_dfs_per_el 3D element stiffness and mass matrices
%AUTHOR
%	Paul Wilcox (15-Jul-2024 17:15:15)

%Deal with optional argument about which DOFs to use
if isempty(varargin)
	dofs_to_use = [];
else
	dofs_to_use = varargin{1};
end

%Record the local node numbers of the element stiffness matrices
loc_nd = [1  1  1  1  2  2  2  2];

%Record the local DOFs of the element stiffness matrices
loc_df = [1  2  3  4  1  2  3  4];

%Get the DOFs if not specified
if isempty(dofs_to_use)
	dofs_to_use = unique(loc_df);
end

%If any inputs blank, return at this point with just the loc_nd and loc_df
if isempty(nds) || isempty(els) || isempty(D) || isempty(rho)
	el_K = [];
	el_M = [];
	el_C = [];
	el_Q = [];
	[loc_nd, loc_df] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use);
	return
end

%Temporary matrices of nodal coordinates to save time
nds_1_1 = nds(els(:, 1), 1);
nds_1_2 = nds(els(:, 1), 2);
nds_2_1 = nds(els(:, 2), 1);
nds_2_2 = nds(els(:, 2), 2);

%Jacobian
J = zeros(size(els, 1), 1, 1);
J(:, 1, 1) = ((nds_1_1 - nds_2_1) .^ 2 + (nds_1_2 - nds_2_2) .^ 2) .^ (1 ./ 2);

%Stiffness matrix
el_K = zeros(size(els, 1), 8, 8);

%Damping matrix
el_C = zeros(size(els, 1), 8, 8);
el_C(:, 1, 4) = nds_2_2 ./ 4 - nds_1_2 ./ 4;
el_C(:, 1, 8) = nds_2_2 ./ 4 - nds_1_2 ./ 4;
el_C(:, 2, 4) = nds_1_1 ./ 4 - nds_2_1 ./ 4;
el_C(:, 2, 8) = nds_1_1 ./ 4 - nds_2_1 ./ 4;
el_C(:, 4, 1) = nds_2_2 ./ 4 - nds_1_2 ./ 4;
el_C(:, 4, 2) = nds_1_1 ./ 4 - nds_2_1 ./ 4;
el_C(:, 4, 5) = nds_2_2 ./ 4 - nds_1_2 ./ 4;
el_C(:, 4, 6) = nds_1_1 ./ 4 - nds_2_1 ./ 4;
el_C(:, 5, 4) = nds_2_2 ./ 4 - nds_1_2 ./ 4;
el_C(:, 5, 8) = nds_2_2 ./ 4 - nds_1_2 ./ 4;
el_C(:, 6, 4) = nds_1_1 ./ 4 - nds_2_1 ./ 4;
el_C(:, 6, 8) = nds_1_1 ./ 4 - nds_2_1 ./ 4;
el_C(:, 8, 1) = nds_2_2 ./ 4 - nds_1_2 ./ 4;
el_C(:, 8, 2) = nds_1_1 ./ 4 - nds_2_1 ./ 4;
el_C(:, 8, 5) = nds_2_2 ./ 4 - nds_1_2 ./ 4;
el_C(:, 8, 6) = nds_1_1 ./ 4 - nds_2_1 ./ 4;

%Mass matrix
el_M = zeros(size(els, 1), 8, 8);

%CRemove unwanted DOFs from element matrices
[loc_nd, loc_df, el_K, el_C, el_M] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use, el_K, el_C, el_M);

end
