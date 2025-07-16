function d = fn_signed_dist_to_bdry(pts, bdry_vtcs, varargin)
%SUMMARY
%   Returns signed (positive exterior) shortest distance of point(s) to 
%   boundary surface described by vertices of edges (2D) or triangular 
%   facets (3D). Dimensionality is determined by number of columns in
%   matrix of point coordinates.
%USAGE
%   d = fn_signed_dist_to_bdry(pts, bdry_vtcs, [bdry_fcs, interior_pt])
%AUTHOR
%   Paul Wilcox (2025)
%INPUTS
%   pts - n_pts x n_dim list of query point coordiantes
%   bdry_vtcs - n_nds x n_dim list of boundary vertex coordinates
%   [bdry_fcs - n_fcs x {2 or 3} list of vertex indices for each edge (2D)
%   or triangular facet (3D). This argument is optional for 2D cases; if
%   omitted it is assumed that bdry_nds describe a closed polygon in order 
%   so bdry_fcs = [1, 2; 2, 3; 3, 4; ... ; n_nds - 1, n_nds; n_nds, 1].
%   This argument is required for 3D cases.
%   [interior_pt] - 1 x n_dim coordinates of point inside boundary surface,
%   which is used to determine correct overall sign of d. If not specified,
%   it is assumed that each edge defiend by nodes [A, B] has exterior on RHS of
%   vector AB (2D) or a facet defined by nodes [A, B, C] with surface 
%   normal calculated as cross product of vectors AB and AC is outwards. If 
%   interior_pt is specified, the vertex ordering of edges or facets 
%   doesn't matter (this function will make them self-consistent and the 
%   position of the interior point will be used to determine the overall 
%   sign of the result (interior = negative). Note that if interior point
%   is not specified, care must be taken to order nodes of each facet
%   correctly otherwise erratic behaviour will be obtained.
%OUTPUTS
%   d - n_pts x 1 signed distance of each point to nearest point on 
%   boundary where sign is negative (interior) or positive (exterior).
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
        d = fn_dist_point_to_bdry_2D_v2(pts, bdry_vtcs, bdry_fcs, interior_pt);
    case 3
        if isempty(bdry_fcs)
            error('3D problems require bdry_fcs to be specified')
        end
        if size(bdry_fcs) ~= 3
            error('bdry_fcs must be n_fcs x 3 for 3D problems')
        end
        d = fn_dist_point_to_bdry_3D(pts, bdry_vtcs, bdry_fcs, interior_pt);
end

end