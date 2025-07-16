function d = fn_dist_point_to_bdry_2D_v2(pts, bdry_vtcs, bdry_fcs, interior_pt)
%SUMMARY
%   Returns signed (positive exterior) shortest distance of point(s) to 
%   boundary surface described by vertices of triangular facets
%USAGE
%   d = fn_dist_point_to_bdry_2D_v2(pts, bdry_nds, bdry_fcs, interior_pt)
%AUTHOR
%   Paul Wilcox (2025)
%INPUTS
%   pts - n_pts x 2 list of query point coordiantes
%   bdry_nds - n_nds x 2 list of boundary vertex coordinates
%   bdry_vtcs - [] or n_fcs x 2 list of vertex indices for each edge
%   interior_pt - [] or 1 x 2 coordinates of point inside boundary surface,
%   which is used to determine correct overall sign of d. If empty,
%   it is assumed that edge specified by vertices [A, B] has exterior on 
%   RHS of vector AB. If interior_pt is specified, the vertex ordering 
%   doesn't matter (this function will make them self-consistent and the 
%   position of the interior point will be used to determine the overall 
%   sign of the result (interior = negative)
%OUTPUTS
%   d - n_pts x 1 signed distance of each point to nearest point on 
%   boundary where sign is negative (interior) or positive (exterior).
%NOTES
%   Formulated to be efficient for checking large numbers of points (i.e.
%   n_pts is large) rather than a large number of facets
%--------------------------------------------------------------------------

n_pts = size(pts, 1);
n_nds = size(bdry_vtcs, 1);

if isempty(bdry_fcs)
    bdry_fcs = [1:n_nds; [2:n_nds, 1]]';
end
n_fcs = size(bdry_fcs, 1);

if ~isempty(interior_pt)
    %In this case, no assumptions are made about order of nodes on each
    %edge defining outward surface normal and they are first shuffled so that
    %they are at least consistently inwards (or outwards)
    %After procedure complete, sign is flipped to get interior_pt at
    %negative distance from surface
    if any(accumarray(bdry_fcs(:), 1) > 2)
        error('A vertex is associated with more than two edges')
    end

    %if vertex appears twice, node order for one of the edges it appears in
    %should be reversed
    for k = 1:n_nds
        [i, j] = find(bdry_fcs == k);
        if i(1) == i(2)
            error('Same vertex is used at both ends of an edge')
        end
        if j(1) == j(2)
            bdry_fcs(i(2), :) = fliplr(bdry_fcs(i(2), :));
        end
    end

    %stick the interior point at the end of pts to be tested (it will be
    %removed later
    pts = [pts; interior_pt];
    n_pts = n_pts + 1;
end


%First get the unit normal vector for each face and the internal
%angle of each vertex
fc_normals = zeros(n_fcs, 2);
fc_vertices = reshape(bdry_vtcs(bdry_fcs(:), :), [size(bdry_fcs), 2]);
for v1 = 1:2
    v2 = mod(v1, 2) + 1;
    a12 = reshape(fc_vertices(:, v2, :) - fc_vertices(:, v1, :), [size(fc_vertices, 1), size(fc_vertices, 2)]); %note cannot use squeeze as that causes bdrys with only 1 face to have first dim collapsed too
    if v1 == 1
        fc_normals = [a12(:, 2), -a12(:, 1)];
    end
end
fc_normals = fc_normals ./ sqrt(sum(fc_normals .^ 2, 2));

%Work out vertices and effective normals for each vertex
nd_normals = zeros(n_nds, 2);
for i = 1:n_nds
    [f, n] = find(bdry_fcs == i);
    for j = 1:numel(f)
        nd_normals(i, :) = nd_normals(i, :) + fc_normals(f(j), :) * 0.5;
    end
end
nd_normals = nd_normals ./ sqrt(sum(nd_normals .^ 2, 2));

%fn_debug_plot(bdry_fcs, bdry_nds, fc_normals, eds, ed_normals, nd_normals)

%Now look in turn for the nearest vertex, edge and face to each point and
%take the one that gives the smallest absolute result as the answer. Sign
%of distance is obtained by sign of dot-product from nearest point with
%effective normal direction.

d = ones(n_pts, 1) * inf;

%Vertices
nds = bdry_vtcs(unique(bdry_fcs(:)), :);
for i = 1:n_nds
    vec = pts - nds(i, :);
    dps = sign(sum(vec .* nd_normals(i,:), 2));
    dps(dps == 0) = 1; %Force sign to be +/1 1, never zero
    r_nds = fn_dist_point_to_point(pts, nds(i, :)) .* dps;
    d = min(d, r_nds, 'ComparisonMethod', 'abs');
end

%Faces
for i = 1:n_fcs
    [r_fcs, alpha, above] = fn_dist_point_to_line(pts, ...
        bdry_vtcs(bdry_fcs(i, 1), :), ...
        bdry_vtcs(bdry_fcs(i, 2), :));

    r_fcs(~above) = inf;
    vec = pts - (bdry_vtcs(bdry_fcs(i, 1), :) + ...
        (bdry_vtcs(bdry_fcs(i, 2), :) - bdry_vtcs(bdry_fcs(i, 1), :)) .* alpha);
    dps = sign(sum(vec .* fc_normals(i,:), 2));
    dps(dps == 0) = 1; %Force sign to be +/1 1, never zero
    r_fcs = r_fcs .* dps;
    d = min(d, r_fcs, 'ComparisonMethod', 'abs');
end

if ~isempty(interior_pt)
    if d(end) > 0
        d = -d;
    end
    d = d(1:end - 1);
end

end
%------------------
%debugging plotting functions

function fn_debug_plot(bdry_fcs, bdry_nds, fc_normals, eds, ed_normals, nd_normals)
arrow_len = sqrt(sum((max(bdry_nds) - min(bdry_nds)) .^ 2)) / 10;
figure;
patch('Faces', bdry_fcs, 'Vertices', bdry_nds,'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
view(3); axis equal; hold on;

for i = 1:size(bdry_fcs, 1)
    fc_cent = mean(bdry_nds(bdry_fcs(i, :), :));
    fn_plot_vec(fc_cent, fc_normals(i,:) * arrow_len, 'r');
end
for i = 1:size(eds, 1)
    ed_cent = mean(bdry_nds(eds(i, :), :));
    fn_plot_vec(ed_cent, ed_normals(i,:) * arrow_len, 'k');
end
for i = 1:size(bdry_nds, 1)
    fn_plot_vec(bdry_nds(i,:), nd_normals(i,:) * arrow_len, 'b');
end
end

function fn_plot_vec(cent, vec, col)
plot3([cent(1), cent(1) + vec(1)], ...
    [cent(2), cent(2) + vec(2)], ...
    [cent(3), cent(3) + vec(3)], ...
    [col, '-']);
plot3(cent(1), cent(2), cent(3), [col, 'o']);
end
