function mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness)
%SUMMARY
%   Adds an absorbing boundary by increasing the element absorbing index
%   proportional to their distance from the specified boundary divided by
%   the specified absorbing boundary layer thickness (i.e. so it reaches
%   one when the distance is equal to the absorbing boundary layer 
%   thickness. The boundary defines the start of the absorbing layer;
%   within the boundary the absorbing index is set to zero.

el_ctrs = fn_calc_element_centres(mod.nds, mod.els);
d = fn_dist_point_to_bdry_2D(el_ctrs, abs_bdry_pts);
in = fn_elements_in_region(mod, abs_bdry_pts);

mod.el_abs_i = d / abs_bdry_thickness;
mod.el_abs_i(mod.el_abs_i < 0) = 0;
mod.el_abs_i(mod.el_abs_i > 1) = 1;
mod.el_abs_i(in) = 0;

end