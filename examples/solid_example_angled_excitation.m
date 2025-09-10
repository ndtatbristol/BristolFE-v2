clearvars -except scripts_to_run
close all;
restoredefaultpath;
addpath(genpath('../code'));

%This example shows how to apply the input excitation as angled forces on
%one edge of the model. Can be applied either parallel (shear) or normal
%(normal) to the surface according to which is selected below:

src_dir = 'shear'; 
%src_dir = 'normal'; 

show_geom_only = 0; %Set to 1 to just show geometry without running model

%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

%Material properties

%Steel
steel_matl_i = 1;
matls{steel_matl_i}.rho = 8900; %Density
matls{steel_matl_i}.D = fn_isotropic_stiffness_matrix(210e9, 0.3); 
matls{steel_matl_i}.col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls{steel_matl_i}.name = 'Steel';

%Element type to use
el_typ_solid = 'CPE3'; 


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
mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);
mod.el_mat_i(:) = steel_matl_i;
el_types = {el_typ_solid};
mod.el_typ_i(:) = find(strcmp(el_types, el_typ_solid));

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
frc_nds = fn_find_nodes_nearest_to_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
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
if ~exist('scripts_to_run') && show_geom_only %suppress graphics when running all scripts for testing
    figure; 
    display_options.draw_elements = 0;
    display_options.node_sets_to_plot(1).nd = steps{1}.load.frc_nds;
    display_options.node_sets_to_plot(1).col = 'r.';
    h_patch = fn_show_geometry(mod, matls, display_options);
    drawnow
    return
end

%--------------------------------------------------------------------------
%RUN THE MODEL

res = fn_FE_entry_point(mod, matls, el_types, steps, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

if ~exist('scripts_to_run') %suppress graphics when running all scripts for testing
    %Show the history output as a function of time - here we sum over all 
    %the nodes/DoFs where displacments were recorded with same weighting as 
    %applied forces to effectively simulate transducer in pulse-echo mode
    figure;
    plot(steps{1}.load.time, steps{1}.load.wts.' * res{1}.dsps);
    xlabel('Time (s)')
    
    if ~isinf(fe_options.field_output_every_n_frames)
        %Animate result
        figure;
        display_options.draw_elements = 0;
        h_patch = fn_show_geometry(mod, matls, display_options);
        anim_options.repeat_n_times = 1;
        fn_run_animation(h_patch, res{1}.fld, anim_options);
    end
end
