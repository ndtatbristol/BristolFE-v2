function mod = fn_add_crack(mod, crack_vtcs, crack_fcs, cod)
%SUMMARY
%   Adds a crack into a 2D or 3D model by identifying nearest element
%   edges/faces and 'splitting' model along them, by duplicating nodes.
%   Default is a zero width crack unless optional Crack Opening 
%   Displacement (COD) is specified in which case the nodes are displaced
%   away from plane of crack
%INPUTS
%   mod - structured variable describing model, containing nodal 
%   coordinates, mod.nds, and element nodes, mod.els. Note that ndim = 
%   size(mod.nds, 2)
%   crack_vtcs - ndim x n_vtcs matrix of coordinates describing vertices of
%   surface that will define crack
%   crack_fcs - n_faces x 2 or 3 matrix of vertex indices for each face
%   (i.a. a line facet in 2D, or triangular facet in 3D). In 2D this
%   parameter can be empty, in which case crack_vtcs are assumed to be
%   consecutive points along a single crack
%   [cod - crack opening displacement, default = 0]
%OUTPUT
%   mod - model with modified nodes and elements.
%--------------------------------------------------------------------------

ndim = size(mod.nds, 2);

if ndim == 3
    [crack_fcs, crack_eds, ~] = fn_consistent_facet_nodes(crack_fcs);
end

%Get signed distance of each element from crack
el_cents = fn_calc_element_centres(mod.nds, mod.els);
d = fn_signed_dist_to_bdry(el_cents, crack_vtcs, crack_fcs);

%identify elements on either side of crack and within tol of crack
[~, max_el_size] = fn_get_min_max_element_sizes(mod);

e0 = abs(d) < max_el_size * 2; %this is purely for computational efficiency to minimuise number of elements considered
ep = (sign(d)) >= 0 & e0;
en = (sign(d)) <  0 & e0;

%Identify the nodes on the crack
[interface_fcs, ~, ~] = fn_find_interface(mod, ep, en);
crack_nd_i = unique(interface_fcs(:));
crack_nds = mod.nds(crack_nd_i, :);

%Get details about each crack_nd relative to crack surface
[~, ~, norm_vecs, type_of_nearest_entity, nearest_entity, bdry_edges] = fn_signed_dist_to_bdry(crack_nds, crack_vtcs, crack_fcs);

switch ndim
    case 2
        %Identify end vertices (ones that only appear once in list of edges)
        tmp = accumarray(bdry_edges(:), 1);
        end_vtcs = find(tmp == 1);
        %Nodes to drop are those that are closest to end vertices rather than
        %other vertices or crack facet
        nodes_to_drop = type_of_nearest_entity == 1 & ismember(nearest_entity, end_vtcs);
    case 3
        %Identify edge of crack surface by edges that only appear once in original
        %list
        [tmp, ~, ic] = unique(sort(crack_eds, 2), 'rows');
        edge_edges = tmp(accumarray(ic, 1) == 1, :);
        edge_edge_i = find(ismember(sort(bdry_edges, 2), edge_edges, 'rows'));

        %Identify vectices on edge of crack (any on edge edges)
        edge_vtc_i = unique(tmp(:));

        %Nodes to drop are those that are closest to edge edges or edge 
        %vertices rather than other vertices, edges or crack facets
        nodes_to_drop = (type_of_nearest_entity == 1 & ismember(nearest_entity, edge_vtc_i)) | ...
            (type_of_nearest_entity == 2 & ismember(nearest_entity, edge_edge_i));

end

crack_nd_i(nodes_to_drop) = [];
crack_nds(nodes_to_drop, :) = [];
norm_vecs(nodes_to_drop, :) = [];

%Apply the crack opening displacement
mod.nds(crack_nd_i, :) = mod.nds(crack_nd_i, :) + norm_vecs .* cod;
crack_nds = crack_nds - norm_vecs .* cod;

%Duplicate crack nodes
new_node_indices = (1:numel(crack_nd_i))' + size(mod.nds,1);
mod.nds = [mod.nds; crack_nds];

%Finally loop through crack_nds and for any occurences in elements on -ve 
%side, change to equivalent new nd
tmp = mod.els(en, :);
for i = 1:numel(crack_nd_i)
    tmp(tmp == crack_nd_i(i)) = new_node_indices(i);
end
mod.els(en, :) = tmp;

end

