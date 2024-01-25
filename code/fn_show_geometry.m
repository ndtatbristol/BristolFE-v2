function h_patch = fn_show_geometry(mod, matls, options)
%SUMMARY
%   Plots the model geometry and returns handle to patches used for each
%   element that can be used as an input argument for subsequent field 
%   animations if desired.

% addpath(genpath('..'));
default_options.offset = [0, 0];
default_options.scale = 1;
default_options.norm_value = []; %empty to normalise to max

options = fn_set_default_fields(options, default_options);

options.matl_cols = zeros(numel(matls), 3);
for i = 1:numel(matls)
    options.matl_cols(i, :) = matls(i).col;
end

if isfield(mod, 'el_mat_i')
    options.el_mat_i = mod.el_mat_i;
else
    options.el_mat_i = ones(size(mod.els, 1), 1);
end


if isfield(mod, 'el_abs_i')
    options.el_abs_i = mod.el_abs_i;
else
    options.el_abs_i = zeros(size(mod.els, 1), 1);
end

h_patch = fn_display_result_v2(mod.nds * options.scale + options.offset, mod.els, options);
end

