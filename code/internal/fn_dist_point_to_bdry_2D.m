function d = fn_dist_point_to_bdry_2D(pts, bdry)
% bdry = [bdry; bdry(1,:)]; %close the loop
%returns shortest distance between point(s) and polygon described by bdry
d = zeros(size(pts, 1), 1);
nb = size(bdry, 1);
for i = 1:size(pts, 1)
    %Find nearest bdry nd
    r0 = sqrt(sum((pts(i, :) - bdry) .^ 2, 2));
    %Get distance to adjacent edges
    [r0, j] = min(r0);
    j1 = rem(j + 1 - 1 + nb, nb) + 1;
    j2 = rem(j - 1 - 1 + nb, nb) + 1;
    [r1, s] = fn_dist_pt_to_line(pts(i, :), bdry(j,:), bdry(j1,:));
    if s < 0 | s > 1
        r1 = inf;
    end
    [r2, s] = fn_dist_pt_to_line(pts(i, :), bdry(j,:), bdry(j2,:));
    if s < 0 | s > 1
        r2 = inf;
    end
    d(i) = min([r0, r1, r2]);
end

end

function [r, s] = fn_dist_pt_to_line(pt, p1, p2)
n = p2 - p1;
d = sqrt(sum(n.^ 2, 2));
n = n / d;
s = sum((pt - p1) .* n, 2);
r = pt - p1 - s * n;
s = s / d;
r = sqrt(sum(r .^ 2, 2));
end