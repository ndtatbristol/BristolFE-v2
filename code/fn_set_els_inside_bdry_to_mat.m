function mod = fn_set_els_inside_bdry_to_mat(mod, bdry_pts, mat_i)
%SUMMARY
%   Sets the material of all elements within boundary to the specified
%   material index

%Work out which elements are in water and which in steel
in = fn_elements_in_region(mod, bdry_pts);

%Change material of those in region to mat
mod.el_mat_i(in) = mat_i;
end