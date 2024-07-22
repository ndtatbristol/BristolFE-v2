clearvars -except scripts_to_run
close all;
restoredefaultpath;
addpath(genpath('../code'));
addpath(genpath('../subdoms'));

%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

abs_bdry_thickness = 1e-3;

steel_matl_i = 1;
%Material properties
main.matls(steel_matl_i).rho = 8900; %Density
main.matls(steel_matl_i).D = fn_isotropic_stiffness_matrix(210e9, 0.3);
main.matls(steel_matl_i).col = hsv2rgb([2/3,0,0.80]); %Colour for display
main.matls(steel_matl_i).name = 'Steel';
main.matls(steel_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

water_matl_i = 2;
main.matls(water_matl_i).rho = 1000;
%For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
%calculated here from ultrasonic velocity (1500) and density (1000)
main.matls(water_matl_i).D = 1500 ^ 2 * 1000;
main.matls(water_matl_i).col = hsv2rgb([0.6,0.5,0.8]);
main.matls(water_matl_i).name = 'Water';
main.matls(water_matl_i).el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid

%Define shape of model
model_size = 10e-3;
bdry_pts = [
    0, 0
    model_size, 0
    model_size, model_size
    0, model_size];

%Define region that will be water
water_bdry_pts = [
    0, 0
    model_size, 0
    model_size, 0.4 * model_size
    0, 0.6 * model_size];

%Define start of absorbing boundary region and its thickness
abs_bdry_pts = [
    abs_bdry_thickness, abs_bdry_thickness
    model_size - abs_bdry_thickness, abs_bdry_thickness
    model_size - abs_bdry_thickness, model_size - abs_bdry_thickness
    abs_bdry_thickness, model_size - abs_bdry_thickness];

%Define a line along which sources will be placed to excite waves
src_end_pts = [0.3, 0.1; 0.7, 0.1] * model_size;
src_dir = 4; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
centre_freq = 5e6;
no_cycles = 4;
max_time = 20e-6;
fe_options.time_pts = 8000;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 12;
safety_factor = 3;

%For animations, set the following to a non-infinite value
fe_options.field_output_every_n_frames = 40;
%--------------------------------------------------------------------------
%PREPARE THE MESH

%Work out element size and Create the nodes and elements of the mesh
el_size = fn_get_suitable_el_size(main.matls, centre_freq, els_per_wavelength);
main.mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%Timestep
main.mod.max_safe_time_step = fn_get_suitable_time_step(main.matls, el_size, safety_factor);
main.mod.design_centre_freq = centre_freq;


%First set material of all elements to steel then set elements inside water
%boundary material to water
main.mod.el_mat_i(:) = steel_matl_i;
main.mod = fn_set_els_inside_bdry_to_mat(main.mod, water_bdry_pts, water_matl_i);

%Add interface elements - this is crucial otherwise there will be no
%coupling between fluid and solid
main.mod = fn_add_fluid_solid_interface_els(main.mod, main.matls);

%Define the absorbing layer
main.mod = fn_add_absorbing_layer(main.mod, abs_bdry_pts, abs_bdry_thickness);

%Define transducers
[main.trans{1}.nds, s] = fn_find_nodes_on_line(main.mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
main.trans{1}.dfs = ones(size(main.trans{1}.nds)) * 4;

%Create a subdomain in the middle with a hole in surface as scatterr
subdomain_size = model_size / 10;
scatterer_size = model_size / 10 * 0.8;
inner_bdry = [-1,-1;-1,1;1,1;1,-1] / 2 * subdomain_size + [1, 1] * model_size / 2;
scat_pts = [-1,0;0,1;1,0;0,-1] / 2 * scatterer_size + [1, 1] * model_size / 2;
main.doms{1}.mod = fn_create_subdomain(main.mod, main.matls, inner_bdry, abs_bdry_thickness);
main.doms{1}.mod = fn_add_scatterer(main.doms{1}.mod, main.matls, scat_pts, water_matl_i);

%Show the mesh
if ~exist('scripts_to_run') %suppress graphics when running all scripts for testing
    figure;
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = main.trans{1}.nds;
    display_options.node_sets_to_plot(1).col = 'r.';
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
        fn_run_animation(h_patches, main.doms{1}.val.trans{1}.fld, anim_options);
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