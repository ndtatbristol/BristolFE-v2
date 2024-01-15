function [animation_data, norm_value] = fn_prepare_animation(mod, fld, options)
% default_options.show_what = 'KE';
default_options.offset = [0, 0];
default_options.scale = 1;
default_options.norm_value = []; %empty to normalise to max
options = fn_set_default_fields(options, default_options);

animation_data.h_patch = fn_display_result_v2(mod.nds * options.scale + options.offset, mod.els, options);

animation_data.v = fld;
if isempty(options.norm_value)
    norm_value = max(animation_data.v, [], 'all');
else
    norm_value = options.norm_value;
end

animation_data.v = animation_data.v / norm_value;

animation_data.max_ti = size(animation_data.v, 2);
animation_data.min_ti = 1;
end

