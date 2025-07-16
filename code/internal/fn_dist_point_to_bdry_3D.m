function d = fn_dist_point_to_bdry_3D(pts, bdry_nds, bdry_fcs, interior_pt)
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
%   it is assumed that each facet is specified with consistent
%   counter-clockwise vertex ordering so that surface normal calculated on
%   cross product of AB and AC is always outwards. If interior_pt is
%   specified, the vertex ordering doesn't matter (this function will make
%   them self-consistent and the position of the interior point will be
%   used to determine the overall sign of the result (interior = negative)
%OUTPUTS
%   d - n_pts x 1 signed distance of each point to nearest point on 
%   boundary where sign is negative (interior) or positive (exterior).
%NOTES
%   Formulated to be efficient for checking large numbers of points (i.e.
%   n_pts is large) rather than a large number of facets
%--------------------------------------------------------------------------

n_pts = size(pts, 1);
n_fcs = size(bdry_fcs, 1);
n_nds = size(bdry_nds, 1);

all_eds = [bdry_fcs(:,1), bdry_fcs(:,2)
    bdry_fcs(:,2), bdry_fcs(:,3)
    bdry_fcs(:,3), bdry_fcs(:,1)];
all_ed_fcs = repmat((1:n_fcs)', [3, 1]);

if ~isempty(interior_pt)
    %In this case, no assumptions are made about order of nodes on each facet
    %defining outward surface normal and they are first shuffled so that
    %they are at least consistently inwards (or outwards)
    %After procedure complete, sign is flipped to get interior_pt at
    %negative distance from surface

    %if edge appears twice, node order for one of the faces it appears in
    %should be reversed
    [~,ia,ic] = unique(sort(all_eds, 2), 'rows');
    occs = accumarray(ic,1);
    if any(occs) > 2
        error('More than two facets at same edge')
    end
    occs = find(occs > 1);
    for i = 1:numel(occs)
        j = find(ic == occs(i));
        if all(all_eds(j(1), :) == all_eds(j(2), :))
            %Flip order of nodes in associated face
            flip_face = all_ed_fcs(j(1));
            bdry_fcs(flip_face,:) = fliplr(bdry_fcs(flip_face,:));
            all_eds(j(1), :) = fliplr(all_eds(j(1), :));
        end
    end
    pts = [pts; interior_pt];
    n_pts = n_pts + 1;
end


%First get the unit normal vector for each face and the internal
%angle of each vertex
fc_normals = zeros(n_fcs, 3);
fc_vertex_weights = zeros(n_fcs, 3);
fc_vertices = reshape(bdry_nds(bdry_fcs(:), :), [size(bdry_fcs), 3]);
for v1 = 1:3
    v2 = mod(v1    , 3) + 1;
    v3 = mod(v1 + 1, 3) + 1;
    a21 = squeeze(fc_vertices(:, v2, :) - fc_vertices(:, v1, :));
    a31 = squeeze(fc_vertices(:, v3, :) - fc_vertices(:, v1, :));
    fc_vertex_weights(:, v1) = real(acos(sum(a21 .* a31, 2) ./ sqrt(sum(a21 .^ 2, 2) .* sum(a31 .^ 2, 2))));
    if v1 == 1
        fc_normals = cross(a21, a31, 2);
    end
end
fc_normals = fc_normals ./ sqrt(sum(fc_normals .^ 2, 2));

%Work out edges and effective normals for each edge
all_eds = sort(all_eds, 2);

[eds, ia, ic] = unique(all_eds, 'rows');
n_eds = size(eds, 1);
ed_normals = zeros(n_eds, 3);
for i = 1:size(all_eds, 1)
    ed_normals(ic(i), :) = ed_normals(ic(i), :) + fc_normals(all_ed_fcs(i), :);
end
ed_normals = ed_normals ./ sqrt(sum(ed_normals .^ 2, 2));

%Work out vertices and effective normals for each vertex
nd_normals = zeros(n_nds, 3);
for i = 1:n_nds
    [f, n] = find(bdry_fcs == i);
    for j = 1:numel(f)
        nd_normals(i, :) = nd_normals(i, :) + fc_normals(f(j), :) * fc_vertex_weights(f(j),n(j));
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
nds = bdry_nds(unique(bdry_fcs(:)), :);
for i = 1:n_nds
    vec = pts - nds(i, :);
    dps = sign(sum(vec .* nd_normals(i,:),2));
    dps(dps == 0) = 1; %Force sign to be +/1 1, never zero
    r_nds = fn_dist_point_to_point(pts, nds(i, :)) .* dps;
    d = min(d, r_nds, 'ComparisonMethod', 'abs');
end

%Edges
for i = 1:n_eds
    [r_eds, alpha, above] = fn_dist_point_to_line(pts, ...
        bdry_nds(eds(i, 1), :), ...
        bdry_nds(eds(i, 2), :));
    r_eds(~above) = inf;
    vec = pts - (bdry_nds(eds(i, 1), :) + (bdry_nds(eds(i, 2), :) - bdry_nds(eds(i, 1), :)) .* alpha);
    dps = sign(sum(vec .* ed_normals(i,:),2));
    dps(dps == 0) = 1; %Force sign to be +/1 1, never zero
    r_eds = r_eds .* dps;
    d = min(d, r_eds, 'ComparisonMethod', 'abs');
end

%Faces
for i = 1:n_fcs
    [r_fcs, alpha, beta, above] = fn_dist_point_to_plane(pts, ...
        bdry_nds(bdry_fcs(i, 1), :), ...
        bdry_nds(bdry_fcs(i, 2), :), ...
        bdry_nds(bdry_fcs(i, 3), :));
    r_fcs(~above) = inf;
    vec = pts - (bdry_nds(bdry_fcs(i, 1), :) + ...
        (bdry_nds(bdry_fcs(i, 2), :) - bdry_nds(bdry_fcs(i, 1), :)) .* alpha + ...
        (bdry_nds(bdry_fcs(i, 3), :) - bdry_nds(bdry_fcs(i, 1), :)) .* beta);
    dps = sign(sum(vec .* fc_normals(i,:),2));
    dps(dps == 0) = 1; %Force sign to be +/1 1, never zero
    r_fcs = r_fcs .* dps;
    d = min(d, r_fcs, 'ComparisonMethod', 'abs');
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
