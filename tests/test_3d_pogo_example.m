clear all
% close all;
% restoredefaultpath;
addpath(genpath('../code'));
show_geom_only = 1;
%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

%Material properties
steel_mat_i = 1;
matls(steel_mat_i).rho = 8900; %Density
%3x3 or 6x6 stiffness matrix of material. Here it is isotropic material and
%fn_isotropic_plane_strain_stiffness_matrix(E, v) converts Young's modulus
%and Poisson's ratio into appropriate 3x3 matrix
matls(steel_mat_i).D = fn_isotropic_stiffness_matrix(210e9, 0.3); 
matls(steel_mat_i).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(steel_mat_i).name = 'Steel';
matls(steel_mat_i).el_typ = 'C3D8R'; %C3D8 8 noded brick

gold_mat_i = 2;
matls(gold_mat_i).rho = 19320; %Density
%3x3 or 6x6 stiffness matrix of material. Here it is isotropic material and
%fn_isotropic_plane_strain_stiffness_matrix(E, v) converts Young's modulus
%and Poisson's ratio into appropriate 3x3 matrix
matls(gold_mat_i).D = fn_isotropic_stiffness_matrix(79e9, 0.3); 
matls(gold_mat_i).col = [83, 69, 22] / 100; %Colour for display
matls(gold_mat_i).name = 'Gold';
matls(gold_mat_i).el_typ = 'C3D8R'; %C3D8 8 noded brick

%Define shape of model
model_size = 10e-3;
abs_layer_thickness = 1e-3;
%corner bts
crnr_pts = [
    0, 0, 0
    1, 1, 1] * model_size;

%Define a line along which sources will be placed to excite waves
src_centre = [0.5, 0.5, 1] * model_size;
src_radius = model_size / 10;

src_dir = 3; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

%Details of input signal
centre_freq = 5e6;
no_cycles = 4;
max_time = 10e-6;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 5;

%The default option is field_output_every_n_frames = inf, which means there
%is no field output. Set to a finite value to get a field output.
fe_options.field_output_every_n_frames = inf;
fe_options.solver = 'pogo';

%--------------------------------------------------------------------------
%PREPARE THE MESH

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh
mod = fn_3d_cubic_structured_mesh(crnr_pts, el_size);

%Quick test of absorbing layers
el_ctrs = fn_calc_element_centres(mod.nds, mod.els);
mod.el_abs_i = zeros(size(mod.el_mat_i));
tmp = (abs_layer_thickness - el_ctrs(:,1)) / abs_layer_thickness;
mod.el_abs_i = tmp .* (tmp > 0);

i = el_ctrs(:, 3) < 2e-3;
mod.el_mat_i(i) = gold_mat_i;

%Quick test of making a void
% i = sqrt(sum((el_ctrs - [0.5, 0.5, 0] * model_size) .^ 2, 2)) < 2e-3;
% mod.els(i, :) = [];
% mod.el_mat_i(i) = [];
% mod.el_abs_i(i) = [];
% [mod.nds, mod.els, ~, ~] = fn_remove_unused_nodes(mod.nds, mod.els);

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
steps{1}.load.frc_nds = fn_find_node_at_point(mod.nds, src_centre, el_size);
steps{1}.load.frc_nds = find(...
    abs(mod.nds(:, 3) - src_centre(3)) < el_size / 2 & ...
    sqrt(sum( (mod.nds(:, 1:2) - src_centre(1:2)) .^ 2, 2 )) < src_radius ...
    );


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
if ~exist('scripts_to_run') && show_geom_only %suppress graphics when running all scripts for testing
    figure;
    display_options.transparency = 0.5;
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, display_options);
    return
end
%--------------------------------------------------------------------------
%RUN THE MODEL

% [res, mats] = fn_FE_entry_point(mod, matls, steps, fe_options);
res = fn_FE_entry_point(mod, matls, steps, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Show the history output as a function of time - here we just sum over all 
%the nodes where displacments were recorded
figure;
plot(steps{1}.load.time, sum(res{1}.dsps, 1));
xlabel('Time (s)')

