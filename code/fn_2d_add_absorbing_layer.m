function mod = fn_2d_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness)
%USAGE
%   mod = fn_2d_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness)
%AUTHOR
%   Paul Wilcox (2025)
%SUMMARY
%   Adds an absorbing boundary by increasing element absorbing indices
%   proportional to their distance from the specified boundary divided by
%   the specified absorbing boundary layer thickness (i.e. so it reaches
%   one when the distance is equal to the absorbing boundary layer 
%   thickness. The boundary defines the start of the absorbing layer;
%   within the boundary the absorbing index is set to zero.
%INPUTS
%   mod - structured variable describing model, containing nodal
%   coordinates, mod.nds, and element nodes, mod.els
%   abs_bdry_pts - coordinates of points defining start absorbing boundary.
%   Does not have to be contained fully within model if absorbing
%   boundaries are only required on certain sides but not others
%   abs_bdry_thickness - absorption index for elements outside boundary
%   will be set to their distance from boundary divided by this value (and
%   capped at unity). Typically this is the distance of the abs_bdry_pts
%   from the edge of the model.
%OUTPUTS
%   mod - structured variable describing model, containing nodal
%   coordinates, mod.nds, element nodes, mod.els, and element absorption
%   index, mod.el_abs_i. 
%NOTES
%   If mod.el_abs_i is present in input mod then the contents will be
%   overwritten by this function. Therefore it should only be called once
%   rather than calling it multiple times to for different parts of an 
%   absorbing boundary.
%--------------------------------------------------------------------------


el_ctrs = fn_calc_element_centres(mod.nds, mod.els);
d = fn_dist_point_to_bdry_2D(el_ctrs, abs_bdry_pts);
in = fn_elements_in_region(mod, abs_bdry_pts);

mod.el_abs_i = d / abs_bdry_thickness;
mod.el_abs_i(mod.el_abs_i < 0) = 0;
mod.el_abs_i(mod.el_abs_i > 1) = 1;
mod.el_abs_i(in) = 0;

end