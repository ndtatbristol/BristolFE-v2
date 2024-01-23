%Test of Global to Local model

clear;
close all;
% restoredefaultpath;
restoredefaultpath;
addpath(genpath(fullfile('..', 'code')));

animate_fname = []; %'AFPAC_v2c';
% fname = 'AFPAC_v2';
fname = [];

show_geom_only = 1;

do_no_defect_cases = 1;
do_defect_cases = 0;
do_validation = 0;
do_direct_injection = 0;
els_per_wavelength = 4; %4 for testing, 6+ for real thing

options.field_output_every_n_frames = 20;
options.movie_mode = 1;
scatterer_of_interest = 3;

% options.field_output_every_n_frames = inf;
% options.movie_mode = 0;


%Material properties (SI units used throughout)
water_velocity = 1500;
water_density = 1000;
air_velocity = 340;
air_density = 1;

%centre frequency (used for element size calc)
centre_freq = 5e6;
number_of_cycles = 5;


%main geometry - overall bounds
model_height = 20e-3;
model_length = 40e-3;
h_wall_thick = 12e-3;
v_wall_thick = 10e-3;
r_water_thick = 5e-3;
cladding_thick = 1e-3;
int_radius = 3e-3;
abs_layer_thick = 2e-3;

%transducer
trans_size = 6e-3;
trans_cent = [abs_layer_thick + trans_size / 2, 15e-3];
trans_angd = 20.6; %angle from horisontal: 20.6 for S and 10.8 to L at 45

time_pts = round(els_per_wavelength / 6 * 20000); %put up to 14k for 6 el/lambda

%Materials
matls(1).rho = 8900;
matls(1).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3);
matls(1).col = hsv2rgb([2/3,0,0.80]);
matls(1).name = 'Steel';

matls(2).rho = 7000;
matls(2).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3);
matls(2).col = [0.8672, 0.9375, 0.1875];
matls(2).name = 'Gold';

matls(3).rho = water_density;
matls(3).D = water_velocity ^ 2 * water_density;
matls(3).col = hsv2rgb([0.6,0.5,0.8]);
matls(3).name = 'Water';

matls(4).rho = air_density;
matls(4).D = air_velocity ^ 2 * water_density;
matls(4).col = [1.0, 1.0, 1.0];
matls(4).name = 'Air';

%Define sub-domains
s = 0;
h = trans_cent(2) - h_wall_thick;

%1 At ray entry point
s = s + 1;
subdomain(s).cent = [trans_cent(1) + h * sind(trans_angd), h_wall_thick];
subdomain(s).inner_rad = 2e-3;

% %2 On back wall
s = s + 1;
subdomain(s).cent = [trans_cent(1) + h * sind(trans_angd) + h_wall_thick, 0];
subdomain(s).inner_rad = 2e-3;

%3 On radius
s = s + 1;
subdomain(s).cent = [trans_cent(1) + h * sind(trans_angd) + 2 * h_wall_thick, h_wall_thick];
subdomain(s).inner_rad = 2e-3;

%--------------------------------------------------------------------------
%Build main model
% main = fn_AFPAC_model( ...
%     matls, centre_freq, els_per_wavelength, model_length, model_height, ...
%     h_wall_thick, v_wall_thick, cladding_thick, abs_layer_thick, ...
%     int_radius, r_water_thick, trans_cent, trans_size, trans_angd, subdomain);
main = fn_AFPAC_model_array( ...
    matls, centre_freq, els_per_wavelength, model_length, model_height, ...
    h_wall_thick, v_wall_thick, cladding_thick, abs_layer_thick, ...
    int_radius, r_water_thick, trans_cent, trans_size, trans_angd, subdomain, 2);

display_options = [];
display_options.matl_cols = [
    hsv2rgb([2/3,0,0.80]), %gray
    [0.8672, 0.9375, 0.1875], %gold
    hsv2rgb([0.6,0.5,0.8]), %nice military blue!
    [1.0, 1.0, 1.0]]; %White for air
display_options.draw_elements = 0;
display_options.interface_el_col = 'k';

% if show_geom_only
%     %Show geometry
%     figure;
%     h_patch = fn_show_GL_geometry(main, options);
%     return
% end

%--------------------------------------------------------------------------

%MAIN model run
time_pts = 2500;
main = fn_run_GL_whole_model(main, time_pts, options);


% figure;
% h_patch = fn_show_geometry(main.mod, main.matls, options)
% fn_run_animation_v2(h_patch, main.res.trans{1}.fld, options);

% return

options.doms_to_run = 1;
options.scats_to_run_in = 1;
options.tx_trans = 2;
main = fn_run_GL_sub_model(main, options);

figure;
h_patch = fn_show_geometry(main.doms{1}.scats{1}.mod, main.matls, options)
while 1
fn_run_animation_v2(h_patch, main.doms{1}.scats{1}.trans{2}.fld, options);
end
return
if fname
    save(fname, 'main', "-v7.3");
end

if do_no_defect_cases || do_defect_cases
    if ~do_defect_cases
        options.scats_to_run_in = 2;
    end
    if ~do_no_defect_cases
        options.scats_to_run_in = 1;
    end
    main = fn_run_scatterer_model(main, options);
end



if fname
    save(fname, 'main', "-v7.3");
end

if options.movie_mode
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
    options.doms_to_run = scatterer_of_interest;
    main = fn_run_validation_models(main, options);
    if options.movie_mode
        d = options.doms_to_run;
        s = 1; %s = 1 is always the defect case
        display_options.interface_el_col = 'k';
        display_options.el_mat_i = main.doms{d}.scats{s}.val.mod.el_mat_i;
        display_options.el_abs_i = main.doms{d}.scats{s}.val.mod.el_abs_i;
        figure;
        animation_data = fn_prepare_animation(main.doms{d}.scats{s}.val.mod.nds, main.doms{d}.scats{s}.val.mod.els, main.doms{d}.scats{s}.val.mats.gl_lookup, main.doms{d}.scats{s}.val.res.tx_rx{1}.f_out, display_options);
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
    main = fn_run_direct_injection_models(main, options);
    if options.movie_mode
        d = scatterer_of_interest;
        s = 1; %s = 1 is always the defect case
        display_options.interface_el_col = 'k';
        display_options.el_mat_i = main.mod.el_mat_i;
        display_options.el_abs_i = main.mod.el_abs_i;
        figure;
        animation_data = fn_prepare_animation(main.mod.nds, main.mod.els, main.mats.gl_lookup, main.doms{d}.scats{s}.di.res.tx_rx{1}.f_out, display_options);
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

if ~options.movie_mode
    d = scatterer_of_interest;
    for d = 1:3
    s = 1;
    ymax = 0.1;
    %plot the A-scans for comparison
    t = main.mod.time * 1e6;
    tmax = max(t);
    u_pristine = fn_convolve(main.res.tx_rx{1}.rx, main.mod.desired_input, 2);
    u_total = fn_convolve(main.doms{d}.scats{s}.res.tx_rx{1}.rx, main.mod.desired_input, 2);
    u_scat = u_total - u_pristine;
    u_val = fn_convolve(main.doms{d}.scats{s}.val.res.tx_rx{1}.rx, main.mod.desired_input, 2);
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

