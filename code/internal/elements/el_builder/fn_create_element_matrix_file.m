function fn_create_element_matrix_file(fname, sym_K, sym_C, sym_M, sym_J, loc_nd, loc_df, no_dims, varargin)
fid = fopen(fname, 'wt');

[~, fn_name] = fileparts(fname);

%Header line
fprintf(fid, ['function [el_K, el_C, el_M, loc_nd, loc_df] = ', fn_name, '(nds, els, D, rho, varargin)\n']);
%Comment lines
fprintf(fid, '%%SUMMARY\n');
fprintf(fid, '%%\tThis function was created automatically by fn_create_element_matrix_file\n');
fprintf(fid, '%%\tand contains code to return the stiffness and mass matrices\n');
fprintf(fid, '%%\tfor multiple elements of the same material and type given by the latter\n');
fprintf(fid, ['%%\tpart of the filename, ', fn_name, '.\n']);
fprintf(fid, '%%INPUTS\n');
fprintf(fid, '%%\tnds - n_nds x n_dims matrix of nodal coordinates\n');
fprintf(fid, '%%\tels - n_els x n_nds_per_el matrix of node indices for each elements\n');
fprintf(fid, '%%\tD - ns x ns material stiffness matrix\n');
fprintf(fid, '%%\trho - material density\n');
fprintf(fid, '%%\t[dofs_to_use = [] - optional string listing the DoFs to use, e.g. ''12''. Use [] for all]\n');
fprintf(fid, '%%OUTPUTS\n');
fprintf(fid, '%%\tel_K, el_C, el_M - n_els x n_dfs_per_el x n_dfs_per_el 3D element stiffness and mass matrices\n');
fprintf(fid, '%%AUTHOR\n');
fprintf(fid, ['%%\tPaul Wilcox (', char(datetime), ')\n']);
fprintf(fid, '\n');

%Deal with any extra lines specified in varargin
for i = 1:length(varargin)
    fprintf(fid, [varargin{i}, '\n']);
end

fprintf(fid, '%%Deal with optional argument about which DOFs to use\n');
fprintf(fid, 'if isempty(varargin)\n\tdofs_to_use = [];\nelse\n\tdofs_to_use = varargin{1};\nend\n\n');

fprintf(fid, '%%Record the local node numbers of the element stiffness matrices\n');
fprintf(fid, ['loc_nd = [', num2str(loc_nd'), '];\n\n']);

fprintf(fid, '%%Record the local DOFs of the element stiffness matrices\n');
fprintf(fid, ['loc_df = [', num2str(loc_df'), '];\n\n']);

fprintf(fid, '%%Get the DOFs if not specified\n');
fprintf(fid, 'if isempty(dofs_to_use)\n\tdofs_to_use = unique(loc_df);\nend\n\n');

fprintf(fid, '%%If any inputs blank, return at this point with just the loc_nd and loc_df\n');
fprintf(fid, 'if isempty(nds) || isempty(els) || isempty(D) || isempty(rho)\n\tel_K = [];\n\tel_M = [];\n\tel_C = [];\n\tel_Q = [];\n\t[loc_nd, loc_df] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use);\n\treturn\nend\n\n');

%Shortcut variables for local nds
fprintf(fid, '%%Temporary matrices of nodal coordinates to save time\n');
un_loc_nd = unique(loc_nd);
for i = 1:numel(un_loc_nd)
    for j = 1:no_dims
        fprintf(fid, 'nds_%i_%i = nds(els(:, %i), %i);\n', un_loc_nd(i), j, un_loc_nd(i), j);
    end
end
fprintf(fid, '\n');

%Jacobian determinant
fprintf(fid, '%%Jacobian\n');
fprintf(fid, fn_format_sym_matrix_for_matlab(sym_J, 'J'));

%Stiffness matrix
fprintf(fid, '%%Stiffness matrix\n');
fprintf(fid, fn_format_sym_matrix_for_matlab(sym_K, 'el_K'));

%Damping matrix
fprintf(fid, '%%Damping matrix\n');
fprintf(fid, fn_format_sym_matrix_for_matlab(sym_C, 'el_C'));

%Mass matrix
fprintf(fid, '%%Mass matrix\n');
fprintf(fid, fn_format_sym_matrix_for_matlab(sym_M, 'el_M'));

%Call function to remove unwanted DOFs in all matrices
fprintf(fid, '%%CRemove unwanted DOFs from element matrices\n');
fprintf(fid, '[loc_nd, loc_df, el_K, el_C, el_M] = fn_remove_dofs_from_el_matrices(loc_nd, loc_df, dofs_to_use, el_K, el_C, el_M);\n');

%End line and close
fprintf(fid, '\nend\n');
fclose(fid);


end