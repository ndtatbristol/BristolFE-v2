function h_patch = fn_show_geometry(mod, matls, options)
addpath(genpath('..'));
default_options.offset = [0, 0];
default_options.scale = 1;
default_options.norm_value = []; %empty to normalise to max

options = fn_set_default_fields(options, default_options);

options.matl_cols = zeros(numel(matls), 3);
for i = 1:numel(matls)
    options.matl_cols(i, :) = matls(i).col;
end
options.el_mat_i = mod.el_mat_i;
options.el_abs_i = mod.el_abs_i;
h_patch = fn_display_result_v2(mod.nds * options.scale + options.offset, mod.els, options);
end

