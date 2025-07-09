function [d, alpha, above] = fn_dist_point_to_line(pts, a, c)
%SUMMARY
%   Calculates distance of point(s) in from lines(s) in 2D or 3D.
%USAGE
%   [d, alpha, above] = fn_dist_point_to_line(pts, a, c)
%INPUTS
%   pts - n_pts x n_dims matrix of pts (n_dims = 2 or 3)
%   a and c  - n_planes x n_dims matrices of two points for each line.
%   Note that only one of n_pts and n_planes can be more than one.
%OUTPUTS
%   d - distance of point(s) to line(s)
%   alpha describe the nearest point(s), q, on the line(s) to each 
%   point p:
%       q = c + alpha * (a - c) 
%   above - is true if q is between a and c.
n_lines = size(c, 1);
n_pts = size(pts, 1);
if n_pts > 1 && n_lines > 1
    error('Must be either one point or one line')
end

a2 = dot(a - c, a - c, 2);

d  = zeros(n_pts, n_lines);
alpha = zeros(n_pts, n_lines);
beta  = zeros(n_pts, n_lines);
ac = a - c;
for i = 1: n_lines
    pa = sum((pts - c(i, :)) .* ac(i, :), 2).';
    alpha(:, i) = pa ./ a2(i);
    d(:, i) = sqrt(sum(((pts - c(i, :)) - alpha(:, i) .* ac(i, :)) .^ 2, 2));
end

above = (alpha >= 0) & (alpha <= 1);

end
