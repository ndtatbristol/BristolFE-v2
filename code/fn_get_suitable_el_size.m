function el_size = fn_get_suitable_el_size(matls, nominal_cent_freq, els_per_wavelength)
%SUMMARY
%    This estimates the shortest wavelength likely to be encountered
%    anywhere in the model at the nominal centre frequency which it will
%    run at, and returns a suggested size of element based on this and the
%    specified number of elements per wavelength.

[~, min_vel] = fn_estimate_max_min_vels(matls);
lambda_min = min_vel / nominal_cent_freq;

el_size = lambda_min / els_per_wavelength;
end
