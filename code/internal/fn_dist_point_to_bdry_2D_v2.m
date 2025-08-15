function [d, nearest_pts, norm_vecs, type_of_nearest_entity, nearest_entity, bdry_edges] = fn_dist_point_to_bdry_2D_v2(pts, bndry_nds, bdry_edges, interior_pt)
%SUMMARY
%   Returns signed (positive exterior) shortest distance of point(s) to 
%   boundary surface described by vertices of triangular facets
%USAGE
%   d = fn_dist_point_to_bdry_2D_v2(pts, bdry_nds, bdry_edges, interior_pt)
%   [d, nearest_pts] = fn_dist_point_to_bdry_2D_v2(pts, bdry_nds, bdry_edges, interior_pt)
%   [d, nearest_pts, norm_vecs] = fn_dist_point_to_bdry_2D_v2(pts, bdry_nds, bdry_edges, interior_pt)
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
%   nearest_pts - n_pts x 2 matrix of coordinates of nearest point on
%   boundary associated with each point
%   norm_vecs - n_pts x 2 matrix of unit vectors of boundary surface 
%   normal at each nearest_pt
%   nearest_entity - n_pts x 1 matrix describing what part of boundary is
%   nearest to each point (1 = vertex, 2 = edge)
%NOTES
%   Formulated to be efficient for checking large numbers of points (i.e.
%   n_pts is large) rather than a large number of facets
%--------------------------------------------------------------------------

n_pts = size(pts, 1);
n_nds = size(bndry_nds, 1);
n_dims = 2;

if isempty(bdry_edges)
    bdry_edges = [1:n_nds; [2:n_nds, 1]]';
end
n_fcs = size(bdry_edges, 1);

if ~isempty(interior_pt)
    %In this case, no assumptions are made about order of nodes on each
    %edge defining outward surface normal and they are first shuffled so that
    %they are at least consistently inwards (or outwards)
    %After procedure complete, sign is flipped to get interior_pt at
    %negative distance from surface
    if any(accumarray(bdry_edges(:), 1) > 2)
        error('A vertex is associated with more than two edges')
    end

    %if vertex appears twice, node order for one of the edges it appears in
    %should be reversed
    for k = 1:n_nds
        [i, j] = find(bdry_edges == k);
        if i(1) == i(2)
            error('Same vertex is used at both ends of an edge')
        end
        if j(1) == j(2)
            bdry_edges(i(2), :) = fliplr(bdry_edges(i(2), :));
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
fc_vertices = reshape(bndry_nds(bdry_edges(:), :), [size(bdry_edges), 2]);
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
    [f, n] = find(bdry_edges == i);
    for j = 1:numel(f)
        nd_normals(i, :) = nd_normals(i, :) + fc_normals(f(j), :) * 0.5;
    end
end
nd_normals = nd_normals ./ sqrt(sum(nd_normals .^ 2, 2));

%fn_debug_plot(bdry_edges, bdry_nds, fc_normals, eds, ed_normals, nd_normals)

%Now look in turn for the nearest vertex, edge and face to each point and
%take the one that gives the smallest absolute result as the answer. Sign
%of distance is obtained by sign of dot-product from nearest point with
%effective normal direction.

d = ones(n_pts, 1) * inf;
nearest_pts = zeros(n_pts, n_dims);
norm_vecs = zeros(n_pts, n_dims);
type_of_nearest_entity = zeros(n_pts, 1);
nearest_entity = zeros(n_pts, 1);

%Vertices
nds = bndry_nds(unique(bdry_edges(:)), :);
for i = 1:n_nds
    vec = pts - nds(i, :);
    dps = sign(sum(vec .* nd_normals(i,:), 2));
    dps(dps == 0) = 1; %Force sign to be +/1 1, never zero
    r_nds = fn_dist_point_to_point(pts, nds(i, :)) .* dps;
    
    j = abs(r_nds) < abs(d);
    d(j) = r_nds(j);
    for k = 1:n_dims
        nearest_pts(j, k) = nds(i, k);
        norm_vecs(j, k) = nd_normals(i, k);
    end
    type_of_nearest_entity(j) = 1;
    nearest_entity(j) = i;
end

%Faces (entity = 2)
for i = 1:n_fcs
    [r_fcs, alpha, above] = fn_dist_point_to_line(pts, ...
        bndry_nds(bdry_edges(i, 1), :), ...
        bndry_nds(bdry_edges(i, 2), :));

    r_fcs(~above) = inf;
    nearest_fc_pts = bndry_nds(bdry_edges(i, 2), :) + ...
        (bndry_nds(bdry_edges(i, 1), :) - bndry_nds(bdry_edges(i, 2), :)) .* alpha;
    vec = pts - nearest_fc_pts;
    dps = sign(sum(vec .* fc_normals(i,:), 2));
    dps(dps == 0) = 1; %Force sign to be +/1 1, never zero
    r_fcs = r_fcs .* dps;

    j = abs(r_fcs) < abs(d);
    d(j) = r_fcs(j);
    for k = 1:n_dims
        nearest_pts(j, k) = nearest_fc_pts(j, k);
        norm_vecs(j, k) = fc_normals(i, k);
    end
    type_of_nearest_entity(j) = 2;
    nearest_entity(j) = i;
    
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

function fn_debug_plot(bdry_edges, bdry_nds, fc_normals, eds, ed_normals, nd_normals)
arrow_len = sqrt(sum((max(bdry_nds) - min(bdry_nds)) .^ 2)) / 10;
figure;
patch('Faces', bdry_edges, 'Vertices', bdry_nds,'FaceColor', 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'k');
view(3); axis equal; hold on;

for i = 1:size(bdry_edges, 1)
    fc_cent = mean(bdry_nds(bdry_edges(i, :), :));
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
