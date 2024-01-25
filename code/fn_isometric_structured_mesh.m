function mod = fn_isometric_structured_mesh(bdry_pts, el_size)
%SUMMARY
%   Utility function for generating a isometric structured mesh of triangular
%   elements, that fills the region specified by bdry_nds.
%INPUTS
%   bdry_pts - n_bdry x 2 matrix of coordinates of the n_bdry points that 
%       will define boundary of mesh.
%OUTPUT
%   mod - structured variable containing fields:
%       .nds - n_nds x 2 matrix of coordinates of each of n_nds nodes
%       .els - n_els x 3 matrix of node indices for each of n_els elements
%       .el_mat_i - n_els x 1 matrix of ones as a placeholder for element
%       material indices assigned elsewhere if more than one type of
%       material is used in model
%--------------------------------------------------------------------------
%Figure out a bounding rectangle for whole shape
crnr_pts = [min(bdry_pts); max(bdry_pts)];

%First generate a rectangular mesh of unit height
sin60 = sind(60);

block_size_x = abs(crnr_pts(2, 1) - crnr_pts(1, 1)) + el_size;
block_size_y = (abs(crnr_pts(2, 2) - crnr_pts(1, 2)));

%Work out how many nodes are needed in x and y
nodes_in_x_direction = ceil(block_size_x / el_size) + 1;
nodes_in_y_direction = ceil(block_size_y / el_size / sin60) + 1;

%Work out nodal coordinates
x = linspace(min(crnr_pts(:, 1)) - el_size / 2, max(crnr_pts(:, 1)) + el_size / 2, nodes_in_x_direction);
y = linspace(min(crnr_pts(:, 2)), max(crnr_pts(:, 2)), nodes_in_y_direction);
[node_x_positions, node_y_positions] = meshgrid(x, y);

%Now shuffle rows of x positions back/forward by half an element
node_x_positions(1:2:end, :) = node_x_positions(1:2:end, :) + el_size / 4;
node_x_positions(2:2:end, :) = node_x_positions(2:2:end, :) - el_size / 4;

%Work out node numbers associated with each element (a bit fiddly as you can see) 
node_numbers = reshape([1:numel(node_x_positions)], nodes_in_y_direction, nodes_in_x_direction);

element_node1a = node_numbers(1:2:end-1, 1:end-1);
element_node2a = node_numbers(2:2:end, 2:end);
element_node3a = node_numbers(2:2:end, 1:end-1);

element_node1b = node_numbers(1:2:end-1, 1:end-1);
element_node2b = node_numbers(1:2:end-1, 2:end);
element_node3b = node_numbers(2:2:end, 2:end);

element_node1c = node_numbers(2:2:end-1, 2:end);
element_node2c = node_numbers(3:2:end, 2:end);
element_node3c = node_numbers(3:2:end, 1:end-1);

element_node1d = node_numbers(2:2:end-1, 1:end-1);
element_node2d = node_numbers(2:2:end-1, 2:end);
element_node3d = node_numbers(3:2:end, 1:end-1);

%Final m x 2 matrix of x and y coordinates for each node
mod.nds = [node_x_positions(:), node_y_positions(:)];
%Final n x 3 matrix of 3 node numbers for each element
mod.els = [
    element_node1a(:), element_node2a(:), element_node3a(:)
    element_node1b(:), element_node2b(:), element_node3b(:)
    element_node1c(:), element_node2c(:), element_node3c(:)
    element_node1d(:), element_node2d(:), element_node3d(:)
    ];

%Now remove elements outside original boundary
[in, out] = fn_elements_in_region(mod, bdry_pts);
mod.els(out, :) = [];

%Tidy up by removing unused nodes
[mod.nds, mod.els] = fn_remove_unused_nodes(mod.nds, mod.els);

%Associate each element with a material index = 1
n_els = size(mod.els, 1);
mod.el_mat_i = ones(n_els, 1);
end