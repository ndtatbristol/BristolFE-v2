clear;
close all;
restoredefaultpath;
addpath('../code');

matls(1).rho = 1000;
%For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
%calcualted here from ultrasonic velocity (1500) and density
matls(1).D = 1500 ^ 2 * matls(1).rho;
matls(1).col = hsv2rgb([0.6,0.5,0.8]);
matls(1).name = 'Water';

%Define shape of model
model_size = 10e-3;
bdry_pts = [0, 0; 2, 0; 2, 1; 1, 1] * model_size;

%Define a line along which sources will be placed to excite waves
src_end_pts = [0.3, 0; 0.7, 0] * model_size;
src_dir = 4; %direction of forces applied: 1 = x, 2 = y, 3 = z, 4 = volumetric expansion (for fluids)

%Details of input signal
centre_freq = 5e6;
no_cycles = 4;
max_time = 40e-6;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 12;

%--------------------------------------------------------------------------

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%Associate each element with a material, element type and absorption index
n_els = size(mod.els, 1);
mod.el_mat_i = ones(n_els, 1);
mod.el_abs_i = zeros(n_els, 1);
mod.el_typ_i = repmat({'AC2D3'}, [n_els, 1]);

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
anim_options.repeat_n_times = 1;
h_patch = fn_show_geometry(mod, matls, anim_options);
fn_run_animation_v2(h_patch, res{1}.fld, anim_options);

