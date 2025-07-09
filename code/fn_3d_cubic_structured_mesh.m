function mod = fn_3d_cubic_structured_mesh(crnr_pts, el_size)
%SUMMARY
%   Utility function for generating a 3d structured mesh of cuboidal
%   elements, that fills the cuboidal region specified by crnr_pts.
%INPUTS
%   crnr_pts - 2 x 3 matrix of coordinates of the two corner points that 
%       will define opposite corners of the cuboidal boundary of mesh.
%OUTPUT
%   mod - structured variable containing fields:
%       .nds - n_nds x 3 matrix of coordinates of each of n_nds nodes
%       .els - n_els x 8 matrix of node indices for each of n_els elements
%       .el_mat_i - n_els x 1 matrix of ones as a placeholder for element
%       material indices assigned elsewhere if more than one type of
%       material is used in model
%--------------------------------------------------------------------------
sz = crnr_pts(2,:) - crnr_pts(1,:); %1x3 vector of mesh size in each direction
n = ceil(sz / el_size) + 1; %1x3 vector of number of nodes in each direction
d = sz ./ (n - 1); %1x3  vector of element dims in each direction


i1 = 1:n(1);
i2 = 1:n(2);
i3 = 1:n(3);

[I1, I2, I3] = ndgrid(i1, i2, i3); %3d meshgrid of node indices

%physical coordinate from indices (i = nx3 matrix of indices in 3
%directions; d = 1x3 vector of element size; c = 1x3 coordaintes of mesh 
%corner; returns nx3 of physical coordinates)
fn_xyz = @(i, d, c) (i - 1) .* d + c;
fn_vec = @(m) m(:);

mod.nds = fn_xyz([I1(:), I2(:), I3(:)], d, crnr_pts(1,:));
i = reshape(1:numel(I1), n);
mod.els = [
    fn_vec(i(1:end - 1, 1:end - 1, 1:end - 1)), ... 
    fn_vec(i(2:end,     1:end - 1, 1:end - 1)), ...
    fn_vec(i(2:end,     2:end,     1:end - 1)), ...
    fn_vec(i(1:end - 1, 2:end,     1:end - 1)), ...
    fn_vec(i(1:end - 1, 1:end-1,   2:end)), ... 
    fn_vec(i(2:end,     1:end - 1, 2:end)), ...
    fn_vec(i(2:end,     2:end,     2:end)), ...
    fn_vec(i(1:end - 1, 2:end,     2:end)), ...
    ];

mod.el_mat_i = ones(size(mod.els, 1), 1);
end
