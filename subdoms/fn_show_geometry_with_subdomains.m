function h_patches = fn_show_geometry_with_subdomains(main, display_options)

default_options.sep_frac = 0.1;
default_options.dom_cols = 'rgbcmy';
display_options = fn_set_default_fields(display_options, default_options);


no_doms = numel(main.doms);
max_dom_size = [0, 0];
for d = 1:no_doms
        max_dom_size = max(max_dom_size, max(main.doms{d}.mod.nds) - min(main.doms{d}.mod.nds));
end

main_dims = max(main.mod.nds) - min(main.mod.nds);
sep = sqrt(sum((main_dims .^ 2))) * display_options.sep_frac;
scale = min(main_dims(1) / (no_doms * max_dom_size(1) + (no_doms - 1) * sep), 1);
step = [(max_dom_size(1) + sep), (max_dom_size(2) + sep)] * scale;

%Main
p = 1;
display_options.scale = 1;
display_options.offset = -min(main.mod.nds);
h_patches{p} = fn_show_geometry(main.mod, main.matls, display_options);
display_options.node_sets_to_plot = [];
hold on;
offset = [0, 0];
for d = 1:no_doms
    offset(1) = (d - 1) * step(1);
    p = p + 1;
    offset(2) = -sep;
    display_options.scale = scale;
    display_options.offset = offset - [min(main.doms{d}.mod.nds(:,1)), max(main.doms{d}.mod.nds(:,2))] * scale;
    h_patches{p} = fn_show_geometry(main.doms{d}.mod, main.matls, display_options);
    col = display_options.dom_cols(rem(d - 1, numel(display_options.dom_cols)) + 1);
    tmp = main.doms{d}.mod.inner_bndry_pts;
    tmp = [tmp; tmp(1,:)];
    xy = fn_dom_coord(tmp, display_options.scale, display_options.offset);
    plot(xy(:, 1), xy(:, 2), col);
    xy = fn_main_coord(tmp, -min(main.mod.nds));
    plot(xy(:, 1), xy(:, 2), col);
end

end

function xy = fn_dom_coord(p, scale, offset)
xy = p * scale + offset;
end

function xy = fn_main_coord(p, offset)
xy = p + offset;
end 
