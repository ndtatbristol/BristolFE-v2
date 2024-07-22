clearvars -except scripts_to_run
close all;
restoredefaultpath;
addpath(genpath('../code'));

%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

%This example shows how to apply the input excitation as angled forces on
%one edge of the model

%Material properties
matls(1).rho = 8900; %Density
%3x3 or 6x6 stiffness matrix of material. Here it is isotropic material and
%fn_isotropic_plane_strain_stiffness_matrix(E, v) converts Young's modulus
%and Poisson's ratio into appropriate 3x3 matrix
matls(1).D = fn_isotropic_stiffness_matrix(210e9, 0.3); 
matls(1).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(1).name = 'Steel';
matls(1).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

%Define shape of model
angle_of_top_edge_degs = 20;
model_size = 10e-3;
bdry_pts = [
    0, 0 
    model_size, 0 
    model_size, model_size * (1 - tand(angle_of_top_edge_degs))
    0, model_size];

%Define a line along which sources will be placed to excite waves
src_end_pts = [
    0.3 * model_size, model_size * (1 - 0.3 * tand(angle_of_top_edge_degs))
    0.7 * model_size, model_size * (1 - 0.7 * tand(angle_of_top_edge_degs))];

% src_dir = 'shear'; 
src_dir = 'normal'; 

%Details of input signal
centre_freq = 5e6;
no_cycles = 4;
max_time = 10e-6;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 10;

%The default option is field_output_every_n_frames = inf, which means there
%is no field output. Set to a finite value to get a field output.
fe_options.field_output_every_n_frames = 5;

%--------------------------------------------------------------------------
%PREPARE THE MESH

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh
mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
frc_nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
%for angled forcing, need to force in both DoF at each forcing node, so
%list is twice as long as number of nodes
steps{1}.load.frc_nds = [frc_nds; frc_nds];
steps{1}.load.frc_dfs = [ones(size(frc_nds)); ones(size(frc_nds)) * 2];

%Also provide the time signal for the loading (if this is a vector, it will
%be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
%of different time signals for each frc_nds/frc_dfs. Here the latter is
%needed so that the appropriate components of the desired force can be
%applied to the two DoF at each forcing node. 
time_step = fn_get_suitable_time_step(matls, el_size);
steps{1}.load.time = 0: time_step:  max_time;
steps{1}.load.frcs = fn_gaussian_pulse(steps{1}.load.time, centre_freq, no_cycles);

%Calculate a weighting vector associated with the size of the force 
%component to be applied at each node/DoF. 
switch src_dir
    case 'normal'
        steps{1}.load.wts = ...
            double(steps{1}.load.frc_dfs == 1) * sind(angle_of_top_edge_degs) + ...
            double(steps{1}.load.frc_dfs == 2) * cosd(angle_of_top_edge_degs);
    case 'shear'
        steps{1}.load.wts = ...
            double(steps{1}.load.frc_dfs == 1) * cosd(angle_of_top_edge_degs) + ...
           -double(steps{1}.load.frc_dfs == 2) * sind(angle_of_top_edge_degs);
end

%Also record displacement history at same points (NB there is no reason why
%these have to be same as forcing points)
steps{1}.mon.nds = steps{1}.load.frc_nds;
steps{1}.mon.dfs = steps{1}.load.frc_dfs;

%Show the mesh
if ~exist('scripts_to_run') %suppress graphics when running all scripts for testing
    figure; 
    display_options.draw_elements = 1;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, display_options);
end

%--------------------------------------------------------------------------
%RUN THE MODEL

res = fn_BristolFE_v2(mod, matls, steps, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Show the history output as a function of time - here we sum over all 
%the nodes/DoFs where displacments were recorded with same weighting as 
%applied forces to effectively simulate transducer in pulse-echo mode
figure;
plot(steps{1}.load.time, steps{1}.load.wts.' * res{1}.dsps);
xlabel('Time (s)')

%Animate result
if ~exist('scripts_to_run') %suppress graphics when running all scripts for testing
    figure;
    display_options.draw_elements = 0;
    h_patch = fn_show_geometry(mod, matls, display_options);
    anim_options.repeat_n_times = 1;
    fn_run_animation(h_patch, res{1}.fld, anim_options);
end
