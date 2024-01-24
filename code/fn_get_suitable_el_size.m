function el_size = fn_get_suitable_el_size(matls, nominal_cent_freq, els_per_wavelength)
[~, min_vel] = fn_estimate_max_min_vels(matls);
lambda_min = min_vel / nominal_cent_freq;

el_size = lambda_min / els_per_wavelength;
end
