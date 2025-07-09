function d = fn_dist_point_to_point(a, b)
%SUMMARY
%   Calculates distance from a to b in from a in 2D or 3D.
%USAGE
%   d = fn_dist_point_to_point(a, b)
%INPUTS
%   a and c  - n_pts x n_dims matrices of two points for each line.
%   Either a or c can have a single row.
%OUTPUTS
%   d - distance from a to b

d = sqrt(sum((a - b) .^ 2, 2));
end
