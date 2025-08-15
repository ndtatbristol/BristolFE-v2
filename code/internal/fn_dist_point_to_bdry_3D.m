function [d, nearest_pts, norm_vecs, type_of_nearest_entity, nearest_entity, bdry_edges] = fn_dist_point_to_bdry_3D(pts, bdry_nds, bdry_fcs, interior_pt)
%SUMMARY
%   Returns signed (positive exterior) shortest distance of point(s) to 
%   boundary surface described by vertices of triangular facets
%USAGE
%   d = fn_dist_point_to_bdry_3D(pts, bdry_nds, bdry_fcs, interior_pt)
%AUTHOR
%   Paul Wilcox (2025)
%INPUTS
%   pts - n_pts x 3 list of query point coordiantes
%   bdry_nds - n_nds x 3 list of boundary vertex coordinates
%   bdry_fcs - n_fcs x 3 list of vertex indices for each triangular facet
%   interior_pt - [] or 1 x 3 coordinates of point inside boundary surface,
%   which is used to determine correct overall sign of d. If empty,
%   sign of d will be random.
%OUTPUTS
%   d - n_pts x 1 signed distance of each point to nearest point on 
%   boundary where sign is negative (interior) or positive (exterior).
%   nearest_pts - n_pts x 3 list of nearest points on boundary
%   norm_vecs - n_pts x 3 list of effective surface normal vectors at
%   nearest points
%NOTES
%   Formulated to be efficient for checking large numbers of points (i.e.
%   n_pts is large) rather than a large number of facets
%--------------------------------------------------------------------------

n_pts = size(pts, 1);
n_fcs = size(bdry_fcs, 1);
n_nds = size(bdry_nds, 1);
n_dims = 3;
n_fcs_per_facet = 3;

%Get node ordering for each facet consistent
[bdry_fcs, all_eds, all_ed_fcs] = fn_consistent_facet_nodes(bdry_fcs);

%Stick interior point on end of list of test points if specified (it will
%be removed at end of function)
if ~isempty(interior_pt)
    pts = [pts; interior_pt];
    n_pts = n_pts + 1;
end

%Get the unit normal vector for each face and the internal
%angle of each vertex
fc_normals = zeros(n_fcs, n_dims);
fc_vertex_weights = zeros(n_fcs, n_fcs_per_facet);
fc_vertices = reshape(bdry_nds(bdry_fcs(:), :), [size(bdry_fcs), 3]);
for v1 = 1:3
    v2 = mod(v1    , n_fcs_per_facet) + 1;
    v3 = mod(v1 + 1, n_fcs_per_facet) + 1;
    a21 = reshape(fc_vertices(:, v2, :) - fc_vertices(:, v1, :), [size(fc_vertices, 1), size(fc_vertices, 2)]); %note cannot use squeeze as that causes bdrys with only 1 face to have first dim collapsed too
    a31 = reshape(fc_vertices(:, v3, :) - fc_vertices(:, v1, :), [size(fc_vertices, 1), size(fc_vertices, 2)]); %note cannot use squeeze as that causes bdrys with only 1 face to have first dim collapsed too
    fc_vertex_weights(:, v1) = real(acos(sum(a21 .* a31, 2) ./ sqrt(sum(a21 .^ 2, 2) .* sum(a31 .^ 2, 2))));
    if v1 == 1
        fc_normals = cross(a21, a31, 2);
    end
end
fc_normals = fc_normals ./ sqrt(sum(fc_normals .^ 2, 2));

%Work out edges and effective normals for each edge
all_eds = sort(all_eds, 2);

[bdry_edges, ia, ic] = unique(all_eds, 'rows');
n_eds = size(bdry_edges, 1);
ed_normals = zeros(n_eds, n_dims);
for i = 1:size(all_eds, 1)
    ed_normals(ic(i), :) = ed_normals(ic(i), :) + fc_normals(all_ed_fcs(i), :);
end
ed_normals = ed_normals ./ sqrt(sum(ed_normals .^ 2, 2));

%Work out vertices and effective normals for each vertex
nd_normals = zeros(n_nds, n_dims);
for i = 1:n_nds
    [f, n] = find(bdry_fcs == i);
    for j = 1:numel(f)
        nd_normals(i, :) = nd_normals(i, :) + fc_normals(f(j), :) * fc_vertex_weights(f(j), n(j));
    end
end
nd_normals = nd_normals ./ sqrt(sum(nd_normals .^ 2, 2));

%fn_debug_plot(bdry_fcs, bdry_nds, fc_normals, eds, ed_normals, nd_normals)

%Now look in turn for the nearest vertex, edge and face to each point and
%take the one that gives the smallest absolute result as the answer. Sign
%of distance is obtained by sign of dot-product from nearest point with
%effective normal direction.

d = ones(n_pts, 1) * inf;
nearest_pts = zeros(n_pts, n_dims);
norm_vecs = zeros(n_pts, n_dims);
type_of_nearest_entity = zeros(n_pts, 1);
nearest_entity = zeros(n_pts, 1);

%Vertices (entity = 1)
nds = bdry_nds(unique(bdry_fcs(:)), :);
for i = 1:n_nds
    vec = pts - nds(i, :);
    dps = sign(sum(vec .* nd_normals(i, :), 2));
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

%Edges (entity = 2)
for i = 1:n_eds
    [r_eds, alpha, above] = fn_dist_point_to_line(pts, ...
        bdry_nds(bdry_edges(i, 1), :), ...
        bdry_nds(bdry_edges(i, 2), :));
    r_eds(~above) = inf;
    nearest_ed_pts = (bdry_nds(bdry_edges(i, 1), :) + (bdry_nds(bdry_edges(i, 2), :) - bdry_nds(bdry_edges(i, 1), :)) .* alpha);

    vec = pts - nearest_ed_pts;
    dps = sign(sum(vec .* ed_normals(i,:),2));
    dps(dps == 0) = 1; %Force sign to be +/1 1, never zero
    r_eds = r_eds .* dps;
    
    j = abs(r_eds) < abs(d);
    d(j) = r_eds(j);
    for k = 1:n_dims
        nearest_pts(j, k) = nearest_ed_pts(j, k);
        norm_vecs(j, k) = ed_normals(i, k);
    end
    type_of_nearest_entity(j) = 2;
    nearest_entity(j) = i;
end

%Faces (entity = 3)
for i = 1:n_fcs
    [r_fcs, alpha, beta, above] = fn_dist_point_to_plane(pts, ...
        bdry_nds(bdry_fcs(i, 1), :), ...
        bdry_nds(bdry_fcs(i, 2), :), ...
        bdry_nds(bdry_fcs(i, 3), :));
    r_fcs(~above) = inf;
    nearest_fc_pts = (bdry_nds(bdry_fcs(i, 1), :) + ...
        (bdry_nds(bdry_fcs(i, 2), :) - bdry_nds(bdry_fcs(i, 1), :)) .* alpha + ...
        (bdry_nds(bdry_fcs(i, 3), :) - bdry_nds(bdry_fcs(i, 1), :)) .* beta);
    vec = pts - nearest_fc_pts;
    dps = sign(sum(vec .* fc_normals(i,:),2));
    dps(dps == 0) = 1; %Force sign to be +/1 1, never zero
    r_fcs = r_fcs .* dps;
    j = abs(r_fcs) < abs(d);
    d(j) = r_fcs(j);
    for k = 1:n_dims
        nearest_pts(j, k) = nearest_fc_pts(j, k);
        norm_vecs(j, k) = fc_normals(i, k);
    end
    type_of_nearest_entity(j) = 3;
    nearest_entity(j) = i;
end

if ~isempty(interior_pt)
    if d(end) > 0
        d = -d
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
