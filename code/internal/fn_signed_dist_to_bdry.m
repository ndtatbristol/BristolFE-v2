function [d, nearest_pts, norm_vecs, type_of_nearest_entity, nearest_entity, bndry_edges] = fn_signed_dist_to_bdry(pts, bdry_vtcs, varargin)
%SUMMARY
%   Returns signed (positive exterior) shortest distance of point(s) to 
%   boundary surface described by vertices of edges (2D) or triangular 
%   facets (3D). Dimensionality is determined by number of columns in
%   matrix of point coordinates.
%USAGE
%   d = fn_signed_dist_to_bdry(pts, bdry_vtcs, [bdry_fcs, interior_pt])
%   [d, nearest_pts] = fn_signed_dist_to_bdry(pts, bdry_vtcs, [bdry_fcs, interior_pt])
%AUTHOR
%   Paul Wilcox (2025)
%INPUTS
%   pts - n_pts x n_dim list of query point coordiantes
%   bdry_vtcs - n_nds x n_dim list of boundary vertex coordinates
%   [bdry_fcs - n_fcs x {2 or 3} list of vertex indices for each edge (2D)
%   or triangular facet (3D). This argument is optional for 2D cases; if
%   empty it is assumed that bdry_nds describe a closed polygon in order 
%   so bdry_fcs = [1, 2; 2, 3; 3, 4; ... ; n_nds - 1, n_nds; n_nds, 1].
%   This argument is required for 3D cases. Note that the order of nodes
%   for each facet does not have to be consistent; the function will force
%   them to be consistent (i.e. so adjacent elements have consistent
%   interior/exterior definitions), but the interior / exterior will be 
%   random if an interior point is not specified.
%   [interior_pt] - 1 x n_dim coordinates of point inside boundary surface,
%   which is used to determine correct overall sign of d. If boundary is 2D
%   closed surface, the signed distance will have the correct sign anyway,
%   even if interior_pt is not specified. However, if the boundary is a 2D 
%   open surface or any 3D surface then the sign of the distance will be 
%   random unless interior_pt is specified.
%OUTPUTS
%   d - n_pts x 1 signed distance of each point to nearest point on 
%   boundary where sign is negative (interior) or positive (exterior).
%   nearest_pts - n_pts x n_dims matrix of coordinates of nearest point on
%   boundary associated with each point
%   norm_vecs - n_pts x n_dims matrix of unit vectors of boundary surface 
%   normal at each nearest_pt
%   type_of_nearest_entity - n_pts x 1 matrix of type of entity that each point
%   is nearest to (1 = vertex, 2 = edge, 3 = face)
%   nearest_entity - n_pts x 1 matrix of index of nearest entity to each point
%   bndry_edges - n_edges x 2 matrix of node numbers of edges describing
%   boundary. In 2D, this will be bdry_fcs id specified; in all other cases
%   the list of edges is generated automatically. This is returned so that
%   nearest_entity can be interpreted when the nearest entity is an edge.
%NOTES
%   Formulated to be efficient for checking large numbers of points (i.e.
%   n_pts is large) rather than a large number of edges / facets
%--------------------------------------------------------------------------

%Input checks
n_dims = size(pts, 2);
if n_dims ~= 2 && n_dims ~=3
    error('Number of dimensions must be 2 or 3')
end

if numel(varargin) < 1
    bdry_fcs = [];
else
    bdry_fcs = varargin{1};
end
if numel(varargin) < 2
    interior_pt = [];
else
    interior_pt = varargin{2};
    if ~isempty(interior_pt) && size(interior_pt, 2) ~= n_dims
        error('Dimension inconsistency for interior_pt')
    end
end
if size(bdry_vtcs, 2) ~= n_dims
    error('Dimension inconsistency for bdry_nds')
end

%Final error check and call 2d or 3d function
switch n_dims
    case 2
        if ~isempty(bdry_fcs) && size(bdry_fcs, 2) ~= 2
            error('If specified, bdry_fcs mus be n_fcs x 2 for 2D problems')
        end
            [d, nearest_pts, norm_vecs, type_of_nearest_entity, nearest_entity, bndry_edges] = fn_dist_point_to_bdry_2D_v2(pts, bdry_vtcs, bdry_fcs, interior_pt);
    case 3
        if isempty(bdry_fcs)
            error('3D problems require bdry_fcs to be specified')
        end
        if size(bdry_fcs) ~= 3
            error('bdry_fcs must be n_fcs x 3 for 3D problems')
        end
        [d, nearest_pts, norm_vecs, type_of_nearest_entity, nearest_entity, bndry_edges] = fn_dist_point_to_bdry_3D(pts, bdry_vtcs, bdry_fcs, interior_pt);
end

end