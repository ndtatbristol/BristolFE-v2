clear all;
close all;
restoredefaultpath;
addpath(genpath('../code'));

%This also shows how to export field output animation to an mp4 file for
%inclusion in presentations etc. To do this, anim_options.mp4_out is set to 
%the desired mp4 filename (without the mp4 extension):

% anim_options.mp4_out = 'anim_out2';


%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

steel_matl_i = 1;
%Material properties
matls(steel_matl_i).rho = 8900; %Density
%3x3 or 6x6 stiffness matrix of material. Here it is isotropic material and
%fn_isotropic_plane_strain_stiffness_matrix(E, v) converts Young's modulus
%and Poisson's ratio into appropriate 3x3 matrix
matls(steel_matl_i).D = fn_isotropic_stiffness_matrix(210e9, 0.3); 
matls(steel_matl_i).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(steel_matl_i).name = 'Steel';
matls(steel_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

[C_L, C_S] = fn_estimate_max_min_vels(matls(steel_matl_i));

water_matl_i = 2;
matls(water_matl_i).rho = 1000;
%For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
%calcualted here from ultrasonic velocity (1500) and density (1000)
matls(water_matl_i).D = 1500 ^ 2 * 1000;
matls(water_matl_i).col = hsv2rgb([0.6,0.5,0.8]);
matls(water_matl_i).name = 'Water'; 
matls(water_matl_i).el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid

abs_bdry_thickness = 1e-3;
%Define shape of model
model_size_x = 14e-3;
model_size_z = 10e-3;
bdry_pts = [
    0, -abs_bdry_thickness
    model_size_x, -abs_bdry_thickness
    model_size_x, model_size_z+2e-3 
    0, model_size_z+2e-3];

%Define region that will be water
water_bdry_pts = [
    0, -abs_bdry_thickness 
    model_size_x,  -abs_bdry_thickness
    model_size_x, model_size_z*0.5
    0, model_size_z*0.5];

water_bdry_pts2 = [
    0, model_size_z
    model_size_x,  model_size_z
    model_size_x, model_size_z+2e-3
    0, model_size_z+2e-3];

%Define a line along which sources will be placed to excite waves
src_end_pts = [0.5*model_size_x - 1e-3, 0; 
               0.5*model_size_x + 1e-3, 0];
src_dir = 4; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
centre_freq = 5e6;
no_cycles = 4;
max_time = 15e-6;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 30;

%The default option is field_output_every_n_frames = inf, which means there
%is no field output. Set to a finite value to get a field output.
fe_options.field_output_every_n_frames = inf;

%--------------------------------------------------------------------------
%PREPARE THE MESH

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%First set material of all elements to steel then set elements inside water 
%boundary material to water
mod.el_mat_i(:) = steel_matl_i;
mod = fn_set_els_inside_bdry_to_mat(mod, water_bdry_pts, water_matl_i);
mod = fn_set_els_inside_bdry_to_mat(mod, water_bdry_pts2, water_matl_i);

%Add interface elements - this is crucial otherwise there will be no
%coupling between fluid and solid
mod = fn_add_fluid_solid_interface_els(mod, matls);

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
steps{1}.load.frc_nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
steps{1}.load.frc_dfs = ones(size(steps{1}.load.frc_nds)) * src_dir;

%Also provide the time signal for the loading (if this is a vector, it will
%be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
%of different time signals for each frc_nds/frc_dfs
time_step = fn_get_suitable_time_step(matls, el_size);
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

%Define start of absorbing boundary region and its thickness

abs_bdry_pts = [
    abs_bdry_thickness, 0
    model_size_x - abs_bdry_thickness, 0
    model_size_x - abs_bdry_thickness, model_size_z+1e-3
    abs_bdry_thickness, model_size_z+1e-3];

%Define the absorbing layer
mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);

%Show the mesh
figure; 
display_options.draw_elements = 0;
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
h_patch = fn_show_geometry(mod, matls, display_options);

%--------------------------------------------------------------------------
%RUN THE MODEL

res = fn_BristolFE_v2(mod, matls, steps, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Show the history output as a function of time - here we just sum over all 
%the nodes where displacments were recorded
figure;
plot(steps{1}.load.time, sum(res{1}.dsps));
xlabel('Time (s)')

%Animate result
if ~isinf(fe_options.field_output_every_n_frames)
    figure;
    display_options.draw_elements = 0; %makes it easier to see waves if element edges not drawn
    h_patch = fn_show_geometry(mod, matls, display_options);
    anim_options.repeat_n_times = 1;
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end

time = steps{1}.load.time;
time = time - (no_cycles / 2) / centre_freq;
save('outputs/time2','time')
AScan = sum(res{1}.dsps);
save('outputs/AScan2','AScan')

figure;
plot(time*1e6,AScan)
xlim([0,13])