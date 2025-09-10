function mod = fn_3d_add_crack(mod, crack_vtcs, crack_fcs, varargin)
%USAGE
%   mod = fn_3d_add_crack(mod, crack_vtcs, crack_fcs [, cod])
%SUMMARY
%   Adds a crack into a 3D model by identifying nearest element
%   edges/faces and 'splitting' model along them, by duplicating nodes.
%   Default is a zero width crack unless optional Crack Opening
%   Displacement (COD) is specified in which case the nodes are displaced
%   away from plane of crack
%INPUTS
%   mod - structured variable describing model, containing nodal
%   coordinates, mod.nds, and element nodes, mod.els. Note that ndim =
%   size(mod.nds, 2)
%   crack_vtcs - n_vtcs x 3 matrix of coordinates describing vertices of
%   surface made of triangular facets that will define crack
%   crack_fcs - n_faces x 3 matrix of vertex indices for each triangular
%   facet in 3D
%   [cod - crack opening displacement, default = 0]
%OUTPUT
%   mod - model with modified nodes and elements.
%--------------------------------------------------------------------------
if isempty(varargin)
    cod = 0;
else
    cod = varargin{1};
end

if size(mod.nds, 2) ~= 3
    error('This function is for 3D models. Use fn_2d_add_crack for 2D models')
end

mod = fn_add_crack(mod, crack_vtcs, crack_fcs, cod);

end

