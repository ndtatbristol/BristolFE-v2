function h_patch = fn_show_GL_geometry(main, options)

default_options.sep_frac = 0.1;
default_options.dom_cols = 'rgbcmy';
options = fn_set_default_fields(options, default_options);


no_doms = numel(main.doms);
max_scats = 0;
max_dom_size = [0, 0];
for d = 1:no_doms
    max_scats = max(max_scats, numel(main.doms{d}.scats));
    for s = 1:numel(main.doms{d}.scats)
        max_dom_size = max(max_dom_size, max(main.doms{d}.scats{s}.mod.nds) - min(main.doms{d}.scats{s}.mod.nds));
    end
end

main_dims = max(main.mod.nds) - min(main.mod.nds);
sep = sqrt(sum((main_dims .^ 2))) * options.sep_frac
scale = main_dims(1) / (no_doms * max_dom_size(1) + (no_doms - 1) * sep);
step = [(max_dom_size(1) + sep), (max_dom_size(2) + sep)] * scale;

%Main
p = 1;
options.scale = 1;
options.offset = -min(main.mod.nds);
h_patch{p} = fn_show_geometry(main.mod, main.matls, options);
hold on;
offset = [0, 0];
for d = 1:no_doms
    max_scats = max(max_scats, numel(main.doms{d}.scats));
    offset(1) = (d - 1) * step(1);
    for s = 1:numel(main.doms{d}.scats)
        p = p + 1;
        offset(2) = -(s - 1) * step(2) - sep;
        options.scale = scale;
        options.offset = offset - [min(main.doms{d}.scats{s}.mod.nds(:,1)), max(main.doms{d}.scats{s}.mod.nds(:,2))] * scale;
        h_patch{p} = fn_show_geometry(main.doms{d}.scats{s}.mod, main.matls, options);
        col = options.dom_cols(rem(d - 1, numel(options.dom_cols)) + 1);
        xy = fn_dom_coord(main.doms{d}.mod.bndry_pts, options.scale, options.offset);
        plot(xy(:, 1), xy(:, 2), col);
        xy = fn_main_coord(main.doms{d}.mod.bndry_pts, -min(main.mod.nds));
        plot(xy(:, 1), xy(:, 2), col);
    end
end
end

function xy = fn_dom_coord(p, scale, offset)
xy = p * scale + offset;
end

function xy = fn_main_coord(p, offset)
xy = p + offset;
end 