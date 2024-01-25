function mod = fn_rectangular_structured_mesh(bdry_pts, el_size)
%SUMMARY
%   Utility function for generating a structured mesh of right isocoles triangular
%   elements, that fills the region specified by bdry_nds. However,
%   fn_isometric_structured_mesh is probably superior!
%INPUTS
%   bdry_pts - n_bdry x 2 matrix of coordinates of the n_bdry points that 
%       will define boundary of mesh.
%OUTPUT
%   mod - structured variable containing fields:
%       .nds - n_nds x 2 matrix of coordinates of each of n_nds nodes
%       .els - n_els x 3 matrix of node indices for each of n_els elements

%--------------------------------------------------------------------------
%Figure out a bounding rectangle for whole shape
crnr_pts = [min(bdry_pts); max(bdry_pts)];

block_size_x = abs(crnr_pts(2, 1) - crnr_pts(1, 1));
block_size_y = abs(crnr_pts(2, 2) - crnr_pts(1, 2));


%Work out how many nodes are needed in x and y
nodes_in_x_direction = ceil(block_size_x / el_size) + 1;
nodes_in_y_direction = ceil(block_size_y / el_size) + 1;

%Work out nodal coordinates
x = linspace(min(crnr_pts(:, 1)), max(crnr_pts(:, 1)), nodes_in_x_direction);
y = linspace(min(crnr_pts(:, 2)), max(crnr_pts(:, 2)), nodes_in_y_direction);
[node_x_positions, node_y_positions] = meshgrid(x, y);

%Work out node numbers associated with each element (a bit fiddly as you can see) 
node_numbers = reshape([1:numel(node_x_positions)], nodes_in_y_direction, nodes_in_x_direction);
element_node1a = node_numbers(1:end-1, 1:end-1);
element_node2a = node_numbers(2:end, 2:end);
element_node3a = node_numbers(2:end, 1:end-1);
element_node1b = node_numbers(1:end-1, 1:end-1);
element_node2b = node_numbers(1:end-1, 2:end);
element_node3b = node_numbers(2:end, 2:end);

%Final m x 2 matrix of x and y coordinates for each node
mod.nds = [node_x_positions(:), node_y_positions(:)];
%Final n x 3 matrix of 3 node numbers for each element
mod.els = [[element_node1a(:), element_node2a(:), element_node3a(:)];[element_node1b(:), element_node2b(:), element_node3b(:)]];

%Now remove elements outside original boundary
[in, out] = fn_elements_in_region(mod, bdry_pts);
mod.els(out, :) = [];

%Tidy up by removing unused nodes
[mod.nds, mod.els] = fn_remove_unused_nodes(mod.nds, mod.els);

%Associate each element with a material index = 1 and absorption index = 0 
%to start with
n_els = size(mod.els, 1);
mod.el_mat_i = ones(n_els, 1);
mod.el_abs_i = zeros(n_els, 1);


end