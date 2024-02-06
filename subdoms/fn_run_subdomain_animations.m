function fn_run_subdomain_animations(main, h_patches, anim_options)
default_options.pause_value = 0.001;
default_options.wave_color = [1,1,1];
default_options.mp4_out = [];
default_options.frame_rate = 48;
default_options.db_range = [-40, 0];
default_options.min_ti = 1;
default_options.max_ti = inf;
default_options.ti_step = 1;
default_options.norm_val = [];
default_options.repeat_n_times = 1;
anim_options = fn_set_default_fields(anim_options, default_options);

for t = 1:numel(main.res.trans)
    for d = 1:numel(main.doms)
        fld{1} = main.res.trans{t}.fld;
        if isfield(main.doms{d}, 'res') && ...
                isfield(main.doms{d}.res, 'trans') && ...
                t <= numel(main.doms{d}.res.trans) && ...
                isfield(main.doms{d}.res.trans{t}, 'fld')
            fld{d + 1} = main.doms{d}.res.trans{t}.fld;
        else
            fld{d + 1} = [];
        end
    end
    fn_run_animation(h_patches, fld, anim_options);
end

end