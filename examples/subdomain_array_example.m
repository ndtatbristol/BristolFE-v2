clearvars -except scripts_to_run
close all;
restoredefaultpath;
addpath(genpath('../code'));
addpath(genpath('../subdoms'));

%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

abs_bdry_thickness = 2e-3;

steel_matl_i = 1;
%Material properties
main.matls(steel_matl_i).rho = 8900; %Density
main.matls(steel_matl_i).D = fn_isotropic_stiffness_matrix(210e9, 0.3);
main.matls(steel_matl_i).col = hsv2rgb([2/3,0,0.80]); %Colour for display
main.matls(steel_matl_i).name = 'Steel';
main.matls(steel_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

%Define shape of model
model_size = 20e-3;
bdry_pts = [
    0, 0
    model_size, 0
    model_size, model_size
    0, model_size];

%Define start of absorbing boundary region and its thickness
abs_bdry_pts = [
    abs_bdry_thickness, abs_bdry_thickness
    model_size - abs_bdry_thickness, abs_bdry_thickness
    model_size - abs_bdry_thickness, model_size - abs_bdry_thickness
    abs_bdry_thickness, model_size - abs_bdry_thickness];


%Define array
no_els = 3;
pitch = 0.5e-3;
array_depth = 2e-3;
centre = [model_size / 2, array_depth];

%Details of input signal
centre_freq = 5e6;
no_cycles = 4;
max_time = 20e-6;
fe_options.time_pts = 2000;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 10;
safety_factor = 3;

%The default option is field_output_every_n_frames = inf, which means there
%is no field output. Set to a finite value to get a field output. Note that
%in subdomain models, requesting field output causes the main model to be
%executed twice for each transducer element, once to generate the transfer
%functions and once to generate the field output.
fe_options.field_output_every_n_frames = 10;
%--------------------------------------------------------------------------
%PREPARE THE MESH

%Work out element size and Create the nodes and elements of the mesh
el_size = fn_get_suitable_el_size(main.matls, centre_freq, els_per_wavelength);
main.mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%Timestep
main.mod.max_safe_time_step = fn_get_suitable_time_step(main.matls, el_size, safety_factor);
main.mod.design_centre_freq = centre_freq;


main.mod.el_mat_i(:) = steel_matl_i;

%Define the absorbing layer
main.mod = fn_add_absorbing_layer(main.mod, abs_bdry_pts, abs_bdry_thickness);

%Define array
tc = mean(1:no_els);
for t = 1:no_els
    el_start = centre + [pitch * (t - tc - 0.5), 0];
    el_end =   centre + [pitch * (t - tc + 0.5), 0];
    [main.trans{t}.nds, s] = fn_find_nodes_on_line(main.mod.nds, el_start, el_end, el_size / 2);
    main.trans{t}.dfs = ones(size(main.trans{t}.nds)) * 2; %DF 2 is y direction
end

%Create a subdomain in the middle with a hole in surface as scatterr
a = linspace(0, 2*pi, 361)';
subdomain_size = model_size / 10;
scatterer_size = model_size / 10 * 0.8;
inner_bdry = [cos(a), sin(a)] / 2 * subdomain_size + [1, 1] * model_size / 2;
scat_pts = [cos(a), sin(a)] / 2 * scatterer_size + [1, 1] * model_size / 2;
main.doms{1}.mod = fn_create_subdomain(main.mod, main.matls, inner_bdry, abs_bdry_thickness);
main.doms{1}.mod = fn_add_scatterer(main.doms{1}.mod, main.matls, scat_pts, 0);

%Show the mesh
if ~exist('scripts_to_run') %suppress graphics when running all scripts for testing

    figure;
    display_options.draw_elements = 0;
    col = 'rgbmkyc';
    for t = 1:no_els
        display_options.node_sets_to_plot(t).nd = main.trans{t}.nds;
        display_options.node_sets_to_plot(t).col = [col(rem(t, numel(col)) + 1), '.'];
    end
    h_patch = fn_show_geometry_with_subdomains(main, display_options);
end
%--------------------------------------------------------------------------

%Run main model
main = fn_run_main_model(main, fe_options);

%Run sub-domain model
main = fn_run_subdomain_model(main, fe_options);

%Animate results if requested
if ~isinf(fe_options.field_output_every_n_frames)
    figure;
    anim_options.repeat_n_times = 1;
    anim_options.db_range = [-40, 0];
    anim_options.pause_value = 0.001;
    h_patches = fn_show_geometry_with_subdomains(main, anim_options);
    fn_run_subdomain_animations(main, h_patches, anim_options);
end

%Run validation model
fe_options.validation_mode = 1;
main = fn_run_main_model(main, fe_options);

%Animate validation results if requested
if ~exist('scripts_to_run') %suppress graphics when running all scripts for testing
    if ~isinf(fe_options.field_output_every_n_frames)
        figure;
        anim_options.repeat_n_times = 1;
        anim_options.db_range = [-40, 0];
        anim_options.pause_value = 0.001;
        h_patches = fn_show_geometry(main.doms{1}.val_mod, main.matls, anim_options);
        for t = 1:no_els
            fn_run_animation(h_patches, main.doms{1}.val.trans{t}.fld, anim_options);
        end
    end
end

%View the time domain data and compare wih validation
figure;
i = max(find(abs(main.inp.sig) > max(abs(main.inp.sig)) / 1000));
mv = max(abs(sum(main.doms{1}.res.fmc.time_data(i:end,: ), 2)));
plot(main.doms{1}.res.fmc.time, sum(main.doms{1}.res.fmc.time_data, 2) / mv, 'k');
hold on;
plot(main.doms{1}.val.fmc.time, sum(main.doms{1}.val.fmc.time_data, 2) / mv, 'b');
plot(main.doms{1}.res.fmc.time, (sum(main.doms{1}.res.fmc.time_data, 2) - sum(main.doms{1}.val.fmc.time_data, 2)) / mv, 'r');
ylim([-1,1]);
legend('Sub-domain method', 'Validation', 'Difference');




