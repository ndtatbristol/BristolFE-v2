%Test of subdomain model

clear all;
close all;
% restoredefaultpath;
addpath(genpath(fullfile('..', 'code')));
addpath(genpath(fullfile('..', 'subdoms')));

animate_fname = []; %'AFPAC_v2c';
% fname = 'AFPAC_v2';
fname = [];

show_geom_only = 0;

do_defect_cases = 1;
do_validation = 1;
% do_direct_injection = 0;

els_per_wavelength = 4; %4 for testing, 6+ for real thing
safety_factor = 3;
time_pts = 8000;

fe_options.field_output_every_n_frames = 40;
fe_options.field_output_every_n_frames = inf;
fe_options.dof_to_use = [];%[1,2,4];

subdoms_to_do = {'A','B','C'};
% subdoms_to_do = {'A','C'};

%--------------------------------------------------------------------------
main = fn_create_test_subdomain_model(els_per_wavelength, safety_factor, subdoms_to_do);

% %Test robustness of calcs by scrambling node order in subdomain
% d = 1;
% new_nd_i = randperm(size(main.doms{d}.mod.nds, 1))';
% [main.doms{d}.mod.nds, main.doms{d}.mod.els, old_nd_i] = fn_renumber_nodes(main.doms{d}.mod.nds, main.doms{d}.mod.els, new_nd_i);
% main.doms{d}.mod.main_nd_i = main.doms{d}.mod.main_nd_i(new_nd_i);
% main.doms{d}.mod.bdry_lyrs = main.doms{d}.mod.bdry_lyrs(new_nd_i);

if do_defect_cases
    for d = 1:numel(main.doms)
        a = linspace(0,2*pi,12)';
        cent = mean(main.doms{d}.mod.inner_bndry_pts);
        rmax = mean(sqrt(sum((main.doms{d}.mod.inner_bndry_pts - cent) .^ 2,2)));
        r = (rand(numel(a), 1) + 1) / 2 * rmax * 0.8;
        scat_pts = cent + r .* [cos(a), sin(a)];
        if strcmp(subdoms_to_do{d}, 'B')
            scat_matl = 0;
        else
            scat_matl = 3;
        end
        main.doms{d}.mod = fn_add_scatterer(main.doms{d}.mod, main.matls, scat_pts, scat_matl);
        main.doms{d}.mod.int_el_i = fn_elements_in_region(main.doms{d}.mod, main.doms{d}.mod.inner_bndry_pts);
    end
end

if show_geom_only
    %Show geometry
    figure;
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = main.trans{1}.nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    fn_show_geometry_with_subdomains(main, display_options);
    return
end

%--------------------------------------------------------------------------


%Run main model
main = fn_run_main_model(main, time_pts, fe_options);

%Animate main runs only
% figure;
% h_patch = fn_show_geometry(main.mod, main.matls, fe_options);
% for t = 1:numel(main.res.trans)
%     fn_run_animation(h_patch, main.res.trans{t}.fld, fe_options);
% end
% return

fe_options.doms_to_run = [];
main = fn_run_subdomain_model(main, fe_options);

if ~isinf(fe_options.field_output_every_n_frames)
    %Animate results if requested, and that's the end
    figure;
    anim_options.repeat_n_times = 1;
    anim_options.db_range = [-40, 0];
    anim_options.pause_value = 0.001;
    h_patches = fn_show_geometry_with_subdomains(main, anim_options);
    fn_run_subdomain_animations(main, h_patches, anim_options);
else
    %View the time domain data
    figure;
    for d = 1:numel(main.doms)
        subplot(numel(main.doms), 1, d);
        fmc = fn_extract_FMC_from_subdomain(main, d);
        plot(fmc.time, sum(fmc.time_data,2));
        hold on;
        fmc = fn_extract_FMC_from_main(main);
        plot(fmc.time, sum(fmc.time_data,2));
    end

end


if do_validation
    fe_options.validation_mode = 1;
    main = fn_run_main_model(main, time_pts, fe_options);
    if ~isinf(fe_options.field_output_every_n_frames)
        %Animate results if requested, and that's the end
        figure;
        anim_options.repeat_n_times = 1;
        anim_options.db_range = [-40, 0];
        anim_options.pause_value = 0.001;
        h_patches = fn_show_geometry(main.doms{1}.val_mod, main.matls, anim_options);
        fn_run_animation(h_patches, main.doms{1}.val.trans{1}.fld, anim_options);
    else
    %View the time domain data
        figure;
        for d = 1:numel(main.doms)
            subplot(numel(main.doms), 1, d);
            [fmc, val_fmc] = fn_extract_FMC_from_subdomain(main, d);
            semilogy(fmc.time, abs(sum(fmc.time_data,2)));
            hold on;
            semilogy(val_fmc.time, abs(sum(val_fmc.time_data,2)));
            semilogy(val_fmc.time, abs(sum(fmc.time_data,2) - sum(val_fmc.time_data,2)));
        end
    end
end


return



if fname
    save(fname, 'main', "-v7.3");
end

if do_no_defect_cases || do_defect_cases
    if ~do_defect_cases
        fe_options.scats_to_run_in = 2;
    end
    if ~do_no_defect_cases
        fe_options.scats_to_run_in = 1;
    end
    main = fn_run_scatterer_model(main, fe_options);
end



if fname
    save(fname, 'main', "-v7.3");
end

if fe_options.movie_mode
    %This figure shows animation of main defect free model and the defect
    %subdomains below (with defects)
    figure;
    display_options.doms_to_show = [];
    animation_data = fn_AFPAC_display(main, display_options, 1);
    display_options.gain = master_gain;
    display_options.pause_value = 0.001;
    display_options.mp4_out = [animate_fname, '_main.mp4'];
    % fprintf(('Press any key to start animation\n')); pause;
    fn_run_animation(animation_data, display_options);

    %This figure shows special animation of specified subdomain with and without defect present
    if do_no_defect_cases && do_defect_cases
        figure;
        display_options.doms_to_show = scatterer_of_interest;
        display_options.max_sol_E = animation_data{1}.max_sol_E;
        display_options.max_flu_E = animation_data{1}.max_flu_E;
        animation_data = fn_AFPAC_display(main, display_options, 1);
        display_options.gain = master_gain;
        display_options.pause_value = 0.001;
        display_options.mp4_out = sprintf([animate_fname, '_with_without_defect_d%i.mp4'], display_options.doms_to_show);
        % fprintf(('Press any key to start animation\n')); pause;
        fn_run_animation(animation_data, display_options);
    end
end

%Following just produce validation and DI animations for 3rd defect region
if do_validation
    fe_options.doms_to_run = scatterer_of_interest;
    main = fn_run_validation_models(main, fe_options);
    if fe_options.movie_mode
        d = fe_options.doms_to_run;
        d = 1; %s = 1 is always the defect case
        display_options.interface_el_col = 'k';
        display_options.el_mat_i = main.doms{d}.scats{d}.val.mod.el_mat_i;
        display_options.el_abs_i = main.doms{d}.scats{d}.val.mod.el_abs_i;
        figure;
        animation_data = fn_prepare_animation(main.doms{d}.scats{d}.val.mod.nds, main.doms{d}.scats{d}.val.mod.els, main.doms{d}.scats{d}.val.mats.gl_lookup, main.doms{d}.scats{d}.val.res.tx_rx{1}.f_out, display_options);
        hold on;
        i = main.mod.tx_rx{1}.nds;
        plot(main.mod.nds(i,1), main.mod.nds(i,2), 'r.');
        display_options.gain = master_gain;
        display_options.pause_value = 0.001;
        display_options.mp4_out = sprintf([animate_fname, '_val_d%i.mp4'], d);
        % fprintf(('Press any key to start animation\n')); pause;
        fn_run_animation(animation_data, display_options);
    end
end
if do_direct_injection
    main = fn_run_direct_injection_models(main, fe_options);
    if fe_options.movie_mode
        d = scatterer_of_interest;
        d = 1; %s = 1 is always the defect case
        display_options.interface_el_col = 'k';
        display_options.el_mat_i = main.mod.el_mat_i;
        display_options.el_abs_i = main.mod.el_abs_i;
        figure;
        animation_data = fn_prepare_animation(main.mod.nds, main.mod.els, main.mats.gl_lookup, main.doms{d}.scats{d}.di.res.tx_rx{1}.f_out, display_options);
        hold on;
        i = main.mod.tx_rx{1}.nds;
        plot(main.mod.nds(i,1), main.mod.nds(i,2), 'r.');
        display_options.gain = master_gain;
        display_options.pause_value = 0.001;
        display_options.mp4_out = sprintf([animate_fname, '_DI_d%i.mp4'], d);
        % fprintf(('Press any key to start animation\n')); pause;
        fn_run_animation(animation_data, display_options);
    end
end

if ~fe_options.movie_mode
    d = scatterer_of_interest;
    for d = 1:3
        d = 1;
        ymax = 0.1;
        %plot the A-scans for comparison
        t = main.mod.time * 1e6;
        tmax = max(t);
        u_pristine = fn_convolve(main.res.tx_rx{1}.rx, main.mod.desired_input, 2);
        u_total = fn_convolve(main.doms{d}.scats{d}.res.tx_rx{1}.rx, main.mod.desired_input, 2);
        u_scat = u_total - u_pristine;
        u_val = fn_convolve(main.doms{d}.scats{d}.val.res.tx_rx{1}.rx, main.mod.desired_input, 2);
        mv = max(abs(u_pristine)) * ymax;


        figure;
        subplot(5,1,1);
        plot(t, u_pristine / mv, 'm');
        xlim([0,tmax]); ylim([-1,1]);
        set(gca, 'XTickLabel', [])

        subplot(5,1,2);
        plot(t, u_scat / mv, 'g');
        xlim([0,tmax]); ylim([-1,1]);
        set(gca, 'XTickLabel', [])

        subplot(5,1,3);
        plot(t, u_total / mv, 'r');
        xlim([0,tmax]); ylim([-1,1]);
        set(gca, 'XTickLabel', [])

        subplot(5,1,4);
        plot(t, u_val / mv, 'm');
        xlim([0,tmax]); ylim([-1,1]);
        set(gca, 'XTickLabel', [])

        subplot(5,1,5);
        plot(t, (u_total - u_val) / mv, 'b');
        xlim([0,tmax]); ylim([-1,1]);

        set(gcf, 'Position', [680   540   560   338]);
    end
end

if fname
    save(fname, 'main', "-v7.3");
end

