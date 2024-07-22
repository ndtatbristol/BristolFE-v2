function [K, C, M, gl_lookup] = fn_build_global_matrices_v5(nds, els, el_mat_i, el_abs_i, el_typ_i, matls, fe_options)
%SUMMARY
%   Creates global matrices from mesh definitions
%INPUTS
%   nds - n x nd matrix of nodal coordinates. The row number is the node
%   number; columns are the coordinates of the node.
%   els - m x max_nds_per_el matrix of element nodes. The row number is the element
%   number; columns are the node numbers of the nodes for each element with
%   trailing zeros if element has less nodes than max_nds_per_el
%   el_mat_i - m x 1 vector of element material indices (which refer to
%   materials in 'matls' parameter
%   el_abs_i - m x 1 vector of element absorption level from 0 to 1
%   el_typ_i - m x 1 cell array of element types
%   matls - p x 1 structured variable of materials with fields
%       matls(i).name - string giving name of material
%       matls(i).rho - density of material
%       matls(i).D - 6x6 stiffness matrix of material. 
%   [options- structured parameter of options - see default_options below for explanations]
%OUTPUTS
%   K, C, M - global (d*n) x (d*n) stiffness, damping and mass matrices, where d is no
%   of DOF at each node
%   gl_lookup - n x max(dof) matrix of global row/col indices for given
%   (node, dof) pair
%--------------------------------------------------------------------------

default_options.dof_to_use = []; %Blank uses all of available ones for all elements, subset can be set e.g. as [1,2]
%Following relate to how absorbing regions are created by adding damping
%matrix and reducing stiffness matrix to try and preserve acoustic
%impedance
default_options.damping_power_law = 3;
default_options.max_damping = 3.1415e+07;
default_options.max_stiffness_reduction = 0.01;
default_options.interface_damping_factor = 0;

fe_options = fn_set_default_fields(fe_options, default_options);

if ~isvector(el_mat_i)
    error('Element materials must be a vector');
end
if length(el_mat_i) ~= size(els, 1)
    error('Length of element materials vector must equal number of elements');
end
if ~isvector(matls)
    error('Materials structure must be a vector');
end

fprintf('    Global matrix builder v5 (nodes = %i, elements = %i, ', size(nds, 1), size(els, 1));
t1 = clock;

%find unique element types, max DoF per element, and actual DoFs in use
unique_typs = unique(el_typ_i);
[unique_df, max_el_df] = fn_find_dof_in_use_and_max_dof_per_el(unique_typs, fe_options.dof_to_use);
% df_per_nd = numel(un_df);
max_df = max(unique_df);

%Prepare global matrices
no_nds = size(nds, 1);
total_dof = no_nds * max_df;
total_el_dfs = size(els, 1) * max_el_df ^ 2;
gm_i = zeros(total_el_dfs, 1);
gm_j = zeros(total_el_dfs, 1);
Kvec = zeros(total_el_dfs, 1);
Mvec = zeros(total_el_dfs, 1);
Cvec = zeros(total_el_dfs, 1);



%Loop over unique element types
i1 = 1;
for t = 1:numel(unique_typs)
    fn_el_mats = str2func(['fn_el_', unique_typs{t}]);
    el_i = strcmp(el_typ_i, unique_typs{t});

    %Find unique element matls for this type
    un_mat = unique(el_mat_i(el_i));
    
    %Loop over unique matls for this element type
    for m = 1:numel(un_mat)
        el_i2 = el_mat_i == un_mat(m) & el_i;

        if ~any(el_i2)
            %No elements of this type and material so skip to next material
            continue
        end
        if un_mat(m) > 0 
            D = matls(un_mat(m)).D;
            if isfield(matls(un_mat(m)), 'density') %deal with legacy naming
                rho = matls(un_mat(m)).density;
            else
                rho = matls(un_mat(m)).rho;
            end
        else
            D = 0; %For elements with no material e.g. interface
            rho = 0;
        end

        %Following is legacy - should require matls to have 6x6 stiffness
        %matrix
        if size(D, 1) == 3 && size(D, 2) == 3
            D = [D(:, 1:2), zeros(3, 3), D(:,3)];
            D = [D(1:2, :); zeros(3, 6); D(3, :)];
        end
        

        %Get the element stiffness and mass matrices
        [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_mats(nds, els(el_i2, :), D, rho, fe_options.dof_to_use);

        % %convert loc_df into indices starting at 1 (necessary to avoid
        % %wasting tons of space when dealing with acoustic elements with the
        % %only DOF is 4)
        % [~, loc_df] = ismember(loc_df, un_df);


        %Calculate element damping matrices based on absorbing index of each element
        el_C = el_C + el_M .* el_abs_i(el_i2) .^ fe_options.damping_power_law *  fe_options.max_damping;
        el_K = el_K .* exp(log(fe_options.max_stiffness_reduction) .* el_abs_i(el_i2) .^ (fe_options.damping_power_law + 1));

        %Work out where the element matrices will go in the global matrices
        [loc_nd_i, loc_nd_j] = meshgrid(loc_nd, loc_nd);
        nd_i = reshape(els(el_i2, loc_nd_i), size(el_K));
        nd_j = reshape(els(el_i2, loc_nd_j), size(el_K));
        [df_i, df_j] = meshgrid(loc_df, loc_df);
        df_i = permute(df_i, [3,1,2]);
        df_j = permute(df_j, [3,1,2]);
        df_i = repmat(df_i, [size(el_K, 1), 1, 1]);
        df_j = repmat(df_j, [size(el_K, 1), 1, 1]);
        gi_i = (nd_i - 1) * max_df + df_i;
        gi_j = (nd_j - 1) * max_df + df_j;

        %At this point, we are just building a list of coordinates in
        %global matrices and the associated values to put in them - the
        %actual sparse global matrices are much more efficiently assembled
        %doing it this way
        i2 = i1 + numel(gi_i) - 1;
        gm_i(i1:i2) = gi_i(:);
        gm_j(i1:i2) = gi_j(:);
        Kvec(i1:i2) = el_K(:);
        Mvec(i1:i2) = el_M(:);
        Cvec(i1:i2) = el_C(:);
        i1 = i2 + 1;
    end
end

%Now build actual sparse matrices
K = sparse(gm_i(1:i2), gm_j(1:i2), Kvec(1:i2), total_dof, total_dof);
M = sparse(gm_i(1:i2), gm_j(1:i2), Mvec(1:i2), total_dof, total_dof);
C = sparse(gm_i(1:i2), gm_j(1:i2), Cvec(1:i2), total_dof, total_dof);

%Reduce global matrices by knocking out any zero columns
tmp = repmat([1:no_nds], max_df, 1);
gl_nds = tmp(:);
tmp = repmat([1:max_df], 1, no_nds);
gl_dofs = tmp(:);
[gl_nds, gl_dofs, ~, K, M, C] = fn_reduce_global_matrices(gl_nds, gl_dofs, [], K, M, C);


%Produce global lookup matrix (row = node, col = DOF, content =
%global matrix index associated with node and DOF).
gl_lookup = fn_create_fast_lookup(gl_nds, gl_dofs, no_nds, 0);

%Add interface damping if requested
if fe_options.interface_damping_factor
    C = fn_add_interface_damping(els, el_typ_i, K, C, M, gl_lookup, fe_options.interface_damping_factor);
end

fprintf('DOF = %i) .......... completed in %.2f secs\n', size(K, 1), etime(clock, t1));

end