clearvars -except scripts_to_run
% close all;
restoredefaultpath;
addpath(genpath('../code'));

%This general purpose example shows how to do a 2D solid model with the
%options to:
%   - include a fluid region before the solid region to show how to assign
%   different parts of the model to different materials and add interface
%   elements
%   - add optional absorbing boundary on 3 sides to suppress reflections.

include_fluid_region = 1;
add_absorbing_boundary = 1;

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

%Water
water_matl_i = 2;
matls{water_matl_i}.rho = 1000;
%For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
%calculated here from ultrasonic velocity (1500) and density (1000)
matls{water_matl_i}.D = 1500 ^ 2 * 1000;
matls{water_matl_i}.col = hsv2rgb([0.6,0.5,0.8]);
matls{water_matl_i}.name = 'Water'; 

%Element types to use
el_typ_solid = 'CPE3'; 
el_typ_fluid = 'AC2D3'; 

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
abs_bdry_thickness = 1e-3;
abs_bdry_pts = [
    abs_bdry_thickness, 0
    model_size - abs_bdry_thickness, 0
    model_size - abs_bdry_thickness, model_size - abs_bdry_thickness
    abs_bdry_thickness, model_size - abs_bdry_thickness];

%Define a line along which sources will be placed to excite waves
src_end_pts = [0.3, 0; 0.7, 0] * model_size;

%Details of input signal
centre_freq = 5e6;
no_cycles = 5;
max_time = 50e-6;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 10;

%The default option is field_output_every_n_frames = inf, which means there
%is no field output. Set to a finite value to get a field output.
fe_options.field_output_every_n_frames = 10;

%--------------------------------------------------------------------------
%PREPARE THE MESH

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh
mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);
el_types = {el_typ_solid, el_typ_fluid};

%First set material of all elements to steel ...
mod.el_mat_i(:) = steel_matl_i;
mod.el_typ_i(:) = find(strcmp(el_types, el_typ_solid));

%... then set elements inside water boundary material to water
if include_fluid_region
    els_in_water = fn_elements_in_region(mod, water_bdry_pts);
    mod.el_mat_i(els_in_water) = water_matl_i;
    mod.el_typ_i(els_in_water) = find(strcmp(el_types, el_typ_fluid));
    
    %Add interface elements - this is crucial otherwise there will be no
    %coupling between fluid and solid
    [mod, el_types] = fn_add_fluid_solid_interface_els(mod, el_types);
    
    %Set source direction
    src_dir = 4; %volumetric source if fluid region is used as source will be in fluid
else
    %Set source direction
    src_dir = 2; %forcing in y direction if no fluid is used as source will be on surface of solid
end


%Define the absorbing layer
if add_absorbing_boundary
    mod = fn_2d_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);
end

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
steps{1}.load.frc_nds = fn_find_nodes_nearest_to_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
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
end