clear;
close all;
restoredefaultpath;
addpath('../code');


%Material properties
matls(1).rho = 8900; %Density
%3x3 or 6x6 stiffness matrix of material. Here it is isotropic material and
%fn_isotropic_plane_strain_stiffness_matrix(E, v) converts Young's modulus
%and Poisson's ratio into appropriate 3x3 matrix
matls(1).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3); 
matls(1).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(1).name = 'Steel';
matls(1).el_typ = 'CPE3';

matls(2).rho = 1000;
%For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
%calcualted here from ultrasonic velocity (1500) and density
matls(2).D = 1500 ^ 2 * matls(1).rho;
matls(2).col = hsv2rgb([0.6,0.5,0.8]);
matls(2).name = 'Water';
matls(2).el_typ = 'AC2D3';

%Define shape of model
model_size = 10e-3;
bdry_pts = [0, 0; 1, 0; 1, 1; 0, 1] * model_size;

%Define region that will be water
water_bdry_pts = [0, 0; 1, 0; 1, 0.4; 0, 0.6] * model_size;
water_bdry_pts = [0, 0; 1, 0; 1, 0; 0, 1] * model_size;

%Define a line along which sources will be placed to excite waves
src_end_pts = [0.3, 0; 0.7, 0] * model_size;
src_dir = 4; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
centre_freq = 5e6;
no_cycles = 4;
max_time = 20e-6;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 10;

%--------------------------------------------------------------------------

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%Set elements inside water boundary to water (mat_i = 2)
mod = fn_set_els_inside_bdry_to_mat(mod, water_bdry_pts, 2);

% %Work out which elements are in water and which in steel
% [water_els, steel_els] = fn_elements_in_region(mod, water_bdry_pts);
% 
% 
% 
% %Associate each element with a material, element type and absorption index
% n_els = size(mod.els, 1);
% 
% %Make them all water to start with
% mod.el_mat_i = ones(n_els, 1);
% mod.el_typ_i = repmat({'AC2D3'}, [n_els, 1]);
% 
% %Then change material and element type of those in steel region
% mod.el_mat_i(steel_els) = 2;
% mod.el_typ_i(steel_els) = {'CPE3'};
% 
% mod.el_abs_i = zeros(n_els, 1);

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

%Show the mesh
figure; 
display_options.draw_elements = 1;
display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
display_options.node_sets_to_plot(1).col = 'r.';
h_patch = fn_show_geometry(mod, matls, display_options);

%Actually run the FE model!
fe_options.field_output_every_n_frames = 10;
res = fn_BristolFE_v2(mod, matls, steps, fe_options);

%Show the history output as a function of time - here we just sum over all 
%the nodes where displacments were recorded
figure;
plot(steps{1}.load.time, sum(res{1}.dsps));
xlabel('Time (s)')

%Animate result
figure;
display_options.draw_elements = 0;
h_patch = fn_show_geometry(mod, matls, display_options);
anim_options.repeat_n_times = 1;
anim_options.norm_val = 2.4682e+06;
anim_options.pause_value = 0.1;
fn_run_animation(h_patch, res{1}.fld, anim_options);
