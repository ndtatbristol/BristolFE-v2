function pts = fn_2d_random_walk(npts, step_mean, step_std, angle_mean, angle_std)
%USAGE
%   pts = fn_2d_random_walk(npts, step_mean, step_std, angle_mean, angle_std)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Returns 
%INPUTS
%   npts - number of points in output
%   step_mean - mean step size between points
%   step_std - standard deviation of step size between points
%   angle_mean - mean angle selected at each point to make next step
%   angle_std - standard deviation of angle selected at each point to make
%   next step
%OUTPUT
%   pts - npts x 2 matrix of coordinates on the random wall with first
%   point always (0, 0)
%NOTES
%   Angles are in radians. Distribitions for random steps are normal. The
%   inputs step_mean, step_std, angle_mean, and angle_std can all be 
%   either scalars or vectors of npts
%--------------------------------------------------------------------------
angs = angle_mean(:) + randn(npts, 1) .* angle_std(:);
step_sizes = step_mean(:) + randn(npts, 1) .* step_std(:);
dpts = [cos(angs), sin(angs)] .* step_sizes;
pts = cumsum(dpts);
pts = pts(1:npts, :);
pts(1, :) = 0;
end