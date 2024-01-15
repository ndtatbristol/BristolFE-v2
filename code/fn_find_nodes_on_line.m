function [node_list, s, r] = fn_find_nodes_on_line(nodes, p1, p2, tol)

%SUMMARY
%   Returns list of nodes that lie along line (with specified tolerance)
%   defined by its endpoints
%INPUTS
%   nodes - n x 2 matrix of nodal coordinates
%   p1 - 2 element vector of coordinates of one end of line
%   p2 - 2 element vector of coordinates of other end of line
%   tol - distance specifying tolerance of distance to line for nodes
%OUTPUTS
%   node_list - list of node numbers on line
%   s - parametric description of position of node on line, where 0 is at
%   p1 and 1 is at p2
%   r - distance of each node in list from line

%--------------------------------------------------------------------------
%error checks
if size(nodes, 2) ~= 2
    error('Nodes input must be two column matrix');
end
if length(p1(:)) ~= 2 | length(p2(:)) ~= 2
    error('Points defining line must be 2 element vectors');
end
if ~isscalar(tol)
    error('Tolerance must be scalar value');
end

% a = p1(:);
% b = p2(:);
% n = a - b;
% nn = n / sqrt(sum(n .^ 2));
% p = nodes;
% mat_a = ones(size(p,1),1) * a';
% r = (mat_a - p) - ((mat_a - p) * nn) * nn';
% r = sqrt(sum(r .^ 2, 2));
% s = (mat_a - p) * n / (n' * n);
% node_list = find((r <= tol) & (s >= 0 - tol) &(s <= 1 + tol));
% s = s(node_list);

n = p2 - p1;
d = sqrt(sum(n.^ 2, 2));
n = n / d;
s = sum((nodes - p1) .* n, 2);
r = nodes - p1 - s * n;
s = s / d;
r = sqrt(sum(r .^ 2, 2));

node_list = find((r <= tol) & (s >= 0 - tol) &(s <= 1 + tol));
s = s(node_list);
r = r(node_list);
end