function mod = fn_2d_add_crack(mod, crack_vtcs, varargin)
%USAGE
%   mod = fn_2d_add_crack(mod, crack_vtcs, [crack_fcs, cod])
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Adds a crack into a 2D by identifying nearest element
%   edges/faces and 'splitting' model along them, by duplicating nodes.
%   Default is a zero width crack unless optional Crack Opening
%   Displacement (COD) is specified in which case the nodes are displaced
%   away from plane of crack
%INPUTS
%   mod - structured variable describing model, containing nodal
%   coordinates, mod.nds, and element nodes, mod.els.
%   crack_vtcs - n_vtcs x 2 matrix of coordinates describing vertices of
%   surface that will define crack
%   [crack_fcs - n_faces x 2 matrix of vertex indices for each line facet of
%   crack. This parameter can be empty or omitted, in which case crack_vtcs
%   are assumed to be consecutive points along a single crack]
%   [cod - crack opening displacement, default = 0]
%OUTPUT
%   mod - model with modified nodes and elements.
%--------------------------------------------------------------------------
if numel(varargin) < 1
    crack_fcs = [];
else
    crack_fcs = varargin{1};
end
if numel(varargin) < 2
    cod = 0;
else
    cod = varargin{2};
end

if size(mod.nds, 2) ~= 2
    error('This function is for 2D models. Use fn_3d_add_crack for 3D models')
end

%If crack fcs not specified in 2d then crack is assumed to along
%line of vertices in order listed
if isempty(crack_fcs)
    crack_fcs = [1:size(crack_vtcs, 1) - 1; 2:size(crack_vtcs, 1)]';
end

mod = fn_add_crack(mod, crack_vtcs, crack_fcs, cod);

end

