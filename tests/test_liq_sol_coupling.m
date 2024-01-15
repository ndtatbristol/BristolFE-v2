%Acoustic element test

clear;
close all;
restoredefaultpath;
addpath(genpath(fullfile('..', 'code')));
% addpath(genpath(fullfile(mfilename('fullpath'), '..')));
% return

%Material properties (SI units used throughout)
velocity = 1500;
density = 1000;
bulk_modulus = velocity ^ 2 * density;

%centre frequency (used for element size calc)
centre_freq = 5e6;
number_of_cycles = 4;
els_per_wavelength = 16;

%main geometry
model_corners = [
    -5, -8
    5, -8
    5, 8
    -5, 8] * 1e-3 * 4;

abs_layer_thick = 2e-3;

safety_factor = 3;

%For the processing
time_pts = 4000;

field_output_every_n_frames = inf;

trans_width = 5e-3;
trans_in_water = 1;
trans_pos = 5e-3;

plate_thickness = 6e-3;

%--------------------------------------------------------------------------
%Create root mesh with nodes and materials only
bounding_nds = [min(model_corners); max(model_corners)];

%Material
matls(1).rho = 8900;
matls(1).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3);
matls(1).col = hsv2rgb([2/3,0,0.80]);
matls(1).name = 'Steel';

matls(2).rho = density;
matls(2).D = bulk_modulus;
matls(2).col = hsv2rgb([0.6,0.5,0.8]);
matls(2).name = 'Water';

[max_vel, min_vel] = fn_estimate_max_min_vels(matls);

%Element size and max time step
el_size = min_vel / centre_freq / els_per_wavelength;
max_safe_time_step = el_size / max_vel / safety_factor;

%Background mesh
[mod.nds, mod.els] = fn_isometric_structured_mesh(bounding_nds, el_size);
mod.el_abs_i = zeros(size(mod.els, 1), 1);

%First make them all matl 1
mod.el_typ_i = repmat({'CPE3'}, [size(mod.els, 1), 1]);
mod.el_mat_i = ones(size(mod.els, 1), 1);

%Absorbing layer
abs_bdry = [
    model_corners(1,1) + abs_layer_thick, model_corners(1,2) + abs_layer_thick
    model_corners(2,1) - abs_layer_thick, model_corners(2,2) + abs_layer_thick
    model_corners(3,1) - abs_layer_thick, model_corners(3,2) - abs_layer_thick
    model_corners(4,1) + abs_layer_thick, model_corners(4,2) - abs_layer_thick];

el_centres = fn_calc_element_centres(mod.nds, mod.els);
mod.el_abs_i = fn_dist_point_to_bdry_2D(el_centres, abs_bdry) / abs_layer_thick;
[in, out] = fn_elements_in_region(mod.nds, mod.els, abs_bdry);
mod.el_abs_i(in) = 0;
mod.el_abs_i(mod.el_abs_i > 1) = 1;
mod.el_abs_i(mod.el_abs_i < 0) = 0;

%Water
i = abs(el_centres(:, 2)) > plate_thickness / 2;
input_dof = 4;
mod.el_mat_i(i) = 2;
mod.el_typ_i(i) = {'AC2D3'};

%Interface
mod = fn_add_fluid_solid_interface_els(mod);

%Transducer
trans_angd = 0;

%transmitting transducer
trans_cent = [0, trans_pos];
trans1  = trans_cent - trans_width / 2 * [cosd(trans_angd), sind(trans_angd)];
trans2  = trans_cent + trans_width / 2 * [cosd(trans_angd), sind(trans_angd)];
[trans_nds_tx, s] = fn_find_nodes_on_line(mod.nds, trans1, trans2, el_size / 2);

%receiving transducer
trans_cent = [0, -trans_pos];
trans1  = trans_cent - trans_width / 2 * [cosd(trans_angd), sind(trans_angd)];
trans2  = trans_cent + trans_width / 2 * [cosd(trans_angd), sind(trans_angd)];
[trans_nds_rx, s] = fn_find_nodes_on_line(mod.nds, trans1, trans2, el_size / 2);

%Loading
steps.load.time = [0:time_pts - 1] * max_safe_time_step;
steps.load.frcs = fn_gaussian_pulse(steps.load.time, centre_freq, number_of_cycles);
steps.load.frc_nds = trans_nds_tx;
steps.load.frc_dfs = ones(size(trans_nds_tx)) * 4;

%Monitoring
steps.mon.nds = [trans_nds_tx; trans_nds_rx];
steps.mon.dfs = ones(size(steps.mon.nds)) * 4;
steps.mon.field_output_every_n_frames = field_output_every_n_frames;

% hist_gi = [tx_gi; rx_gi];
pe_nds = [ones(size(trans_nds_tx)); zeros(size(trans_nds_rx))];
pc_nds = [zeros(size(trans_nds_tx)); ones(size(trans_nds_rx))];

%Do it!
options = [];
res = fn_BristolFE_v2(mod, matls, steps, options);

%History output
figure;
pe_sig = sum(res{1}.dsp(find(pe_nds),:))';
pc_sig = sum(res{1}.dsp(find(pc_nds),:))';
plot(steps.load.time, [pe_sig, pc_sig])

%Animate field. Note that the intent is for there NOT to be a separate
%function for preparing animation - it is just done after standard show
%geometry function

return
figure;
options = [];
h_patch = fn_show_geometry(mod, matls, options);
options = [];
fn_run_animation_v2(h_patch, res{1}.fld, options);
