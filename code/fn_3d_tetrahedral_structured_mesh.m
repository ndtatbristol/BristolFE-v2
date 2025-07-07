function mod = fn_3d_tetrahedral_structured_mesh(crnr_pts, el_size)
%SUMMARY
%   Utility function for generating a 3d structured mesh of tetrahedral
%   elements, that fills the cuboidal region specified by crnr_pts.
%INPUTS
%   crnr_pts - 2 x 3 matrix of coordinates of the two corner points that 
%       will define opposite corners of the cuboidal boundary of mesh.
%OUTPUT
%   mod - structured variable containing fields:
%       .nds - n_nds x 3 matrix of coordinates of each of n_nds nodes
%       .els - n_els x 4 matrix of node indices for each of n_els elements
%       .el_mat_i - n_els x 1 matrix of ones as a placeholder for element
%       material indices assigned elsewhere if more than one type of
%       material is used in model
%--------------------------------------------------------------------------
%DUH - regular tetrahedra do not tessalate!!!!!!
dx = el_size;
dy = el_size * sqrt(3) / 2;
dz = el_size * sqrt(6) / 3;

nx = 4;
ny = 4;
nz = 2;

ix = 1:nx;
iy = 1:ny;
iz = 1:nz;

[Ix, Iy, Iz] = ndgrid(ix, iy, iz);

%check if integer is even
fn_even = @(q) 1 - rem(q, 2);
%physical coordinate from indices
fn_xyz = @(i, j, k, dx, dy, dz) [...
    (i + fn_even(j) / 2 + fn_even(k) / 2) * dx, ...
    (j + fn_even(k) / 4) * dy, ...
    k * dz];

mod.nds = fn_xyz(Ix(:), Iy(:), Iz(:), dx, dy, dz);
i = reshape(1:nx*ny*nz, [nx, ny, nz]);
tri1 = i(1:nx-1, 1:2:ny-1, 1:2:nz-1);
tri2 = i(2:nx, 1:2:ny-1, 1:2:nz-1);
mod.els = [
    % tri1(:), tri1(:) + 1, tri1(:) + nx, tri1(:) + nx * ny
    tri1(:) + 1, tri1(:) + nx + 1, tri1(:) + 1, tri1(:) + nx * ny
    ];
end
