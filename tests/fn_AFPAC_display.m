function animation_data = fn_AFPAC_display(main, display_options, prepare_animation)

%Main model
if ~isfield(display_options, 'doms_to_show')
    display_options.doms_to_show = [];
end
display_options.el_mat_i = main.mod.el_mat_i;
display_options.el_abs_i = main.mod.el_abs_i;
display_options.offset = [0, 0];
display_options.scale = 1;
display_options.interface_el_col = 'k';

if isempty(display_options.doms_to_show)
    %This is the usual one
    if prepare_animation
        % animation_data{1} = fn_prepare_animation(main.mod.nds, main.mod.els, main.mats.gl_lookup, main.res.tx_rx{1}.f_out, display_options);
        animation_data{1} = fn_prepare_animation(main.mod.nds, main.mod.els, main.mats.gl_lookup, main.doms{3}.scats{1}.di.res.tx_rx{1}.f_out, display_options);
        
        display_options.max_sol_E = animation_data{1}.max_sol_E;
        display_options.max_flu_E = animation_data{1}.max_flu_E;
    else
        animation_data = [];
        fn_display_result_v2(main.mod.nds * display_options.scale + display_options.offset, main.mod.els, display_options);
    end


    %Transducer
    i = main.mod.tx_rx{1}.nds;
    hold on;
    plot(main.mod.nds(i,1), main.mod.nds(i,2), 'r.');

    %Sub-domains
    s = 1;
    for d = 1:numel(main.doms)
        display_options.offset = [(d - 0) * 14e-3 - 3e-3, -2e-3] - max(main.doms{d}.scats{s}.mod.nds);
        display_options.el_mat_i = main.doms{d}.scats{s}.mod.el_mat_i;
        display_options.el_abs_i = main.doms{d}.scats{s}.mod.el_abs_i;
        if prepare_animation
            animation_data{d + 1} = fn_prepare_animation(main.doms{d}.scats{s}.mod.nds, main.doms{d}.scats{s}.mod.els, main.doms{d}.scats{s}.mats.gl_lookup, main.doms{d}.scats{s}.res.tx_rx{1}.f_out, display_options);
        else
            fn_display_result_v2(main.doms{d}.scats{s}.mod.nds * display_options.scale + display_options.offset, main.doms{d}.scats{s}.mod.els, display_options);
        end
        plot(main.doms{d}.mod.bndry_pts(:,1), main.doms{d}.mod.bndry_pts(:,2), 'g');
        plot(main.doms{d}.mod.bndry_pts(:,1) + display_options.offset(1), main.doms{d}.mod.bndry_pts(:,2) + display_options.offset(2), 'g');
    end
else
    %This is special one to show subdomain with and without scat
    d = display_options.doms_to_show;
    for s = 1:numel(main.doms{d}.scats)
        display_options.offset = [(s - 0) * 14e-3 - 3e-3, 0] - max(main.doms{d}.scats{s}.mod.nds);
        display_options.el_mat_i = main.doms{d}.scats{s}.mod.el_mat_i;
        display_options.el_abs_i = main.doms{d}.scats{s}.mod.el_abs_i;
        if prepare_animation
            animation_data{s} = fn_prepare_animation(main.doms{d}.scats{s}.mod.nds, main.doms{d}.scats{s}.mod.els, main.doms{d}.scats{s}.mats.gl_lookup, main.doms{d}.scats{s}.res.tx_rx{1}.f_out, display_options);
        else
            fn_display_result_v2(main.doms{d}.scats{s}.mod.nds * display_options.scale + display_options.offset, main.doms{d}.scats{s}.mod.els, display_options);
        end
        plot(main.doms{d}.mod.bndry_pts(:,1) + display_options.offset(1), main.doms{d}.mod.bndry_pts(:,2) + display_options.offset(2), 'g');
    end
end
end