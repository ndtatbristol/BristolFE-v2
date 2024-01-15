function node = fn_find_node_at_point(nodes, p, tol)
%SUMMARY
%   Finds nearest node to specified point
%INPUTS
%   nodes - n x 2 matrix of nodal coordinates
%   p - 2 element vector of coordinates of point
%   tol - distance specifying maximum distance of acceptable nodes from
%   point
%OUTPUTS
%   node - number of node nearest to specified point. If no node is within
%   tol of point then zero is returned.

%--------------------------------------------------------------------------

r = sqrt((nodes(:,1) - p(1)) .^ 2 + (nodes(:,2) - p(2)) .^ 2);
[rmin, node] = min(r);
if rmin > tol
    node = 0;
end
end