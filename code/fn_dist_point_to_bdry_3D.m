function d = fn_dist_point_to_bdry_3D(pts, bdry_nds, bdry_fcs, varargin)
%returns shortest distance of point(s) to boundary surface described by n_ndsx3
%matrix of nodal coordinates and n_fcsx3 matrix of nodal indices for n_f
%triangular facets

%returns shortest distance between point(s) and polygon described by bdry
if numel(varargin) < 1
    method = 'loop over boundary';
else
    method = varargin{1};
end
n_pts = size(pts, 1);
d = ones(n_pts, 1) * inf;
switch method
    case 'loop over points'
        for i = 1:n_pts
            %Find closest node
            r_nd = sqrt(sum((pts(i, :) - bdry_nds) .^ 2, 2));
            [r_nd, closest_nd_i] = min(r_nd);
            %Find all faces that involve this node
            closest_fc_i = find(sum(ismember(bdry_fcs, closest_nd_i),2) > 0);
            %Find distance to planes of these faces
            [r_fcs, ~, ~, above] = fn_dist_point_to_plane(pts(i, :), ...
                bdry_nds(bdry_fcs(closest_fc_i, 1), :), ...
                bdry_nds(bdry_fcs(closest_fc_i, 2), :), ...
                bdry_nds(bdry_fcs(closest_fc_i, 3), :));
            r_fcs(~above) = inf;
            %Find all edges that involve this node
            ed_nds = unique(reshape(bdry_fcs(closest_fc_i, :), [],1));
            ed_nds = ed_nds(ed_nds ~= closest_nd_i);
            [r_eds, ~, above] = fn_dist_point_to_line(pts(i, :), ...
                bdry_nds(closest_nd_i, :), ...
                bdry_nds(ed_nds, :));
            r_eds(~above) = inf;
            %Answer for this point is minimum of above
            d(i) = min([r_nd;r_eds(:);r_fcs(:)]);
        end
    case 'loop over boundary'
        %nodes
        nds = bdry_nds(unique(bdry_fcs(:)), :);
        for i = 1:size(nds, 1)
            d = min(d, fn_dist_point_to_point(pts, nds(i, :)));
        end
        %edges
        eds = sort([bdry_fcs(:,1), bdry_fcs(:,2)
            bdry_fcs(:,2), bdry_fcs(:,3)
            bdry_fcs(:,3), bdry_fcs(:,1)], 2);
        eds = unique(eds, 'rows');
        for i = 1:size(eds, 1)
            [r_eds, ~, above] = fn_dist_point_to_line(pts, ...
                bdry_nds(eds(i, 1), :), ...
                bdry_nds(eds(i, 2), :));
            r_eds(~above) = inf;
            d = min(d, r_eds);
        end
        %faces
        for i = 1:size(bdry_fcs, 1)
            [r_fcs, ~, ~, above] = fn_dist_point_to_plane(pts, ...
                bdry_nds(bdry_fcs(i, 1), :), ...
                bdry_nds(bdry_fcs(i, 2), :), ...
                bdry_nds(bdry_fcs(i, 3), :));
            r_fcs(~above) = inf;
            d = min(d, r_fcs);
        end


end
end