function [nodes, elements] = fn_hexagonal_structured_mesh(corner_nodes, element_size)
%SUMMARY
%   Utility function for generating a rectangular structured mesh of equilateral triangular
%   elements arranged in a hexagon fashon
%   Top and bottom surfaces wil be flat; left-right surfaces will wiggle
%INPUTS
%   corner_nodes - 2 x 2 matrix of nodal coordinates - each row is
%   coordinate of diagonally-opposite corners

block_size_x = abs(corner_nodes(2, 1) - corner_nodes(1, 1));
block_size_y = abs(corner_nodes(2, 2) - corner_nodes(1, 2));


%Work out how many nodes are needed in x and y
cos30 = cosd(30);
nodes_in_y_direction = ceil(block_size_y / element_size / cos30) + 1;
y = linspace(min(corner_nodes(:, 2)), max(corner_nodes(:, 2)), nodes_in_y_direction);
element_size = (y(2) - y(1))  / cos30;
nodes_in_x_direction = ceil(block_size_x / element_size) + 1;
x = [1:nodes_in_x_direction] * element_size;
x = x - mean(x) + (corner_nodes(2, 1) + corner_nodes(1, 1)) / 2;

%Work out nodal coordinates
[node_x_positions, node_y_positions] = meshgrid(x, y);
node_x_positions(1:2:end,:) = node_x_positions(1:2:end,:) + element_size / 4;
node_x_positions(2:2:end,:) = node_x_positions(2:2:end,:) - element_size / 4;

%Work out node numbers associated with each element (a bit fiddly as you can see) 
node_numbers = reshape([1:numel(node_x_positions)], nodes_in_y_direction, nodes_in_x_direction);
element_node1a = node_numbers(1:2:end-1, 1:end-1);
element_node2a = node_numbers(2:2:end, 2:end);
element_node3a = node_numbers(2:2:end, 1:end-1);
element_node1b = node_numbers(1:2:end-1, 1:end-1);
element_node2b = node_numbers(1:2:end-1, 2:end);
element_node3b = node_numbers(2:2:end, 2:end);


element_node1c = node_numbers(2:2:end-1, 1:end-1);
element_node2c = node_numbers(2:2:end-1, 2:end);
element_node3c = node_numbers(3:2:end, 1:end-1);
element_node1d = node_numbers(3:2:end, 1:end-1);
element_node2d = node_numbers(2:2:end-1, 2:end);
element_node3d = node_numbers(3:2:end, 2:end);


%Final m x 2 matrix of x and y coordinates for each node
nodes = [node_x_positions(:), node_y_positions(:)];
%Final n x 3 matrix of 3 node numbers for each element
elements = [
    [element_node1a(:), element_node2a(:), element_node3a(:)]
    [element_node1b(:), element_node2b(:), element_node3b(:)]
    [element_node1c(:), element_node2c(:), element_node3c(:)]
    [element_node1d(:), element_node2d(:), element_node3d(:)]
    ];

end