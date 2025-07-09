function [d, alpha, beta, above] = fn_dist_point_to_plane(pts, a, b, c)
%SUMMARY
%   Calculates distance of point(s) in n_pts x 3 matrix pts from plane(s) in 
%   which points a, b, and c lie, where a, b, and c are n_planes x 3 matrices.
%   Only one of n_pts and n_planes can be more than one.
%USAGE
%   [d, alpha, beta, above] = fn_dist_point_to_plane(pts, a, b, c)
%INPUTS
%   pts - n_pts x 3 matrix of pts
%   a, b, and c  - n_planes x 3 matrices of three points in each plane.
%   Note that only one of n_pts and n_planes can be more than one.
%OUTPUTS
%   d - distance of point(s) to plane(s)
%   alpha, beta describe the nearest point(s), q, in the plane(s) to each 
%   point p:
%       q = c + alpha * (a - c) + beta * (b - c) 
%   above - is true if q is inside triangle abc.
%--------------------------------------------------------------------------
n_planes = size(a, 1);
n_pts = size(pts, 1);
if n_pts > 1 && n_planes > 1
    error('Must be either one point or one plane')
end



a2 = dot(a - c, a - c, 2);
b2 = dot(b - c, b - c, 2);
ab = dot(a - c, b - c, 2);



% alpha = zeros(n_planes, n_pts);
% beta  = zeros(n_planes, n_pts);
d  = zeros(n_pts, n_planes);
alpha = zeros(n_pts, n_planes);
beta  = zeros(n_pts, n_planes);
for i = 1: n_planes
    pa = sum((pts - c(i, :)) .* (a(i, :) - c(i, :)), 2).';
    pb = sum((pts - c(i, :)) .* (b(i, :) - c(i, :)), 2).';
    tmp = [a2(i), ab(i); ab(i), b2(i)] \ [pa; pb];
    alpha(:, i) = tmp(1, :) .';
    beta(:, i) = tmp(2, :) .';
    d(:, i) = sqrt(sum(((pts - c(i, :)) - alpha(:, i) .* (a(i, :) - c(i, :)) - beta(:, i) .* (b(i, :) - c(i, :))) .^ 2, 2));
end

above = (alpha >= 0) & (beta >= 0) & ((alpha + beta) <= 1);

end
