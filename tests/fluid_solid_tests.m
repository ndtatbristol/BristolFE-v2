%Fluid-solid transmission, reflection test
clear;
close all;
restoredefaultpath;
addpath(genpath('../code'));

%--------------------------------------------------------------------------
%DEFINE THE PROBLEM

src_angle_degs = [0:2:20] + 180; %0 is normal incidence
% src_angle_degs = 30 + 180;

fluid_vel = 1500;
fluid_rho = 1000;
solid_L_vel = 6000;
solid_S_vel = 3000;
solid_rho = 9000;

refl_angle_degs = -src_angle_degs;
refr_L_angle_degs = -asind(sind(src_angle_degs) * solid_L_vel / fluid_vel);
refr_S_angle_degs = -asind(sind(src_angle_degs) * solid_S_vel / fluid_vel);


show_mesh = 1;
fe_options.interface_damping_factor = 10;

steel_matl_i = 1;
%Material properties
matls(steel_matl_i).rho = 8900; %Density
%3x3 or 6x6 stiffness matrix of material. Here it is isotropic material and
%fn_isotropic_plane_strain_stiffness_matrix(E, v) converts Young's modulus
%and Poisson's ratio into appropriate 3x3 matrix
[lambda, mu] = fn_lame_from_velocities_and_density(solid_L_vel, solid_S_vel, solid_rho);
[youngs_modulus, poissons_ratio] = fn_youngs_from_lame(lambda, mu);
matls(steel_matl_i).D = fn_isotropic_stiffness_matrix(youngs_modulus, poissons_ratio); 
matls(steel_matl_i).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(steel_matl_i).name = 'Solid';
matls(steel_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

water_matl_i = 2;
matls(water_matl_i).rho = fluid_rho;
%For fluids, stiffness 'matrix' D is just the scalar bulk modulus,
%calcualted here from ultrasonic velocity (1500) and density (1000)
matls(water_matl_i).D = fluid_vel ^ 2 * 1000;
matls(water_matl_i).col = hsv2rgb([0.6,0.5,0.8]);
matls(water_matl_i).name = 'Fluid'; 
matls(water_matl_i).el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid


%Define shape of model
model_size = 5e-3;
a = linspace(0, 2 * pi, 361)';
bdry_pts = [cos(a), sin(a)] * model_size;

%Define region that will be water
water_bdry_pts = [
    -1, 0
    1, 0
    1, -1
    -1, -1] * model_size;

%absorbing boundary
abs_bdry_thickness = 1e-3;
abs_bdry_pts = [cos(a), sin(a)] * (model_size - abs_bdry_thickness);

%Define source
src_rad = model_size - 2 * abs_bdry_thickness;
src_dir = 4; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)
src_len = 4e-3;
mon_rad = model_size - 3 * abs_bdry_thickness;

%Details of input signal
centre_freq = 5e6;
no_cycles = 4;
max_time = 8e-6;
safety_factor = 1.5;

%Elements per wavelength (higher = more accurate and higher computational cost)
els_per_wavelength = 10;

%The default option is field_output_every_n_frames = inf, which means there
%is no field output. Set to a finite value to get a field output.
fe_options.field_output_every_n_frames = 20;

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

%Add interface elements - this is crucial otherwise there will be no
%coupling between fluid and solid
mod = fn_add_fluid_solid_interface_els(mod, matls);

%Define the absorbing layer
mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
time_step = fn_get_suitable_time_step(matls, el_size, safety_factor);
t = 0: time_step:  max_time;
in_sig = fn_gaussian_pulse(t, centre_freq, no_cycles);
[~, i] = max(abs(fn_hilbert(in_sig(:))));
t0 = t(i);

time_window = [-1, 1] * no_cycles / centre_freq;

for i = 1:numel(src_angle_degs)
    src_end_pts = src_rad * [sind(src_angle_degs(i)), cosd(src_angle_degs(i))] + [-1;1] * [-cosd(src_angle_degs(i)), sind(src_angle_degs(i))] * src_len / 2;
    steps{i}.load.frc_nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
    steps{i}.load.frc_dfs = ones(size(steps{i}.load.frc_nds)) * src_dir;
    %Also provide the time signal for the loading (if this is a vector, it will
    %be applied at all frc_nds/frc_dfs simultaneously; alternatively it can be a matrix
    %of different time signals for each frc_nds/frc_dfs
    
    steps{i}.load.time = t;
    steps{i}.load.frcs = in_sig;

    steps{i}.mon.nds = [];
    steps{i}.mon.dfs = [];
    steps{i}.mon.case = [];

    %Monitoring points for incident wave in water
    mon_end_pts =  mon_rad * [sind(src_angle_degs(i)), cosd(src_angle_degs(i))] + [-1;1] * [-cosd(src_angle_degs(i)), sind(src_angle_degs(i))] * src_len / 2;
    tmp = fn_find_nodes_on_line(mod.nds, mon_end_pts(1, :), mon_end_pts(2, :), el_size / 2);
    steps{i}.mon.nds = [steps{i}.mon.nds; tmp];
    steps{i}.mon.dfs = [steps{i}.mon.dfs; ones(size(tmp)) * src_dir];
    steps{i}.mon.case = [steps{i}.mon.case; ones(size(tmp)) * 1];
    steps{i}.mon.time_window(1, :) = (src_rad - mon_rad) / fluid_vel + t0 + time_window;

    %Monitoring points for reflected waves
    mon_end_pts = mon_rad * [sind(refl_angle_degs(i)), cosd(refl_angle_degs(i))] + [-1;1] * [-cosd(refl_angle_degs(i)), sind(refl_angle_degs(i))] * src_len / 2;
    tmp = fn_find_nodes_on_line(mod.nds, mon_end_pts(1, :), mon_end_pts(2, :), el_size / 2);
    steps{i}.mon.nds = [steps{i}.mon.nds; tmp];
    steps{i}.mon.dfs = [steps{i}.mon.dfs; ones(size(tmp)) * src_dir];
    steps{i}.mon.case = [steps{i}.mon.case; ones(size(tmp)) * 2];
    steps{i}.mon.time_window(2, :) = (src_rad + mon_rad) / fluid_vel + t0 + time_window;

    %Monitoring points for refracted L waves
    mon_end_pts = mon_rad * [sind(refr_L_angle_degs(i)), cosd(refr_L_angle_degs(i))] + [-1;1] * [-cosd(refr_L_angle_degs(i)), sind(refr_L_angle_degs(i))] * src_len / 2;
    tmp = fn_find_nodes_on_line(mod.nds, mon_end_pts(1, :), mon_end_pts(2, :), el_size / 2);
    steps{i}.mon.nds = [steps{i}.mon.nds; tmp; tmp]; %need both DoF at each node in solid
    steps{i}.mon.dfs = [steps{i}.mon.dfs; ones(size(tmp)); ones(size(tmp)) * 2]; %need both DoF at each node in solid
    steps{i}.mon.case = [steps{i}.mon.case; [ones(size(tmp)); ones(size(tmp))] * 3]; %need both DoF at each node in solid
    steps{i}.mon.time_window(3, :) = (src_rad / fluid_vel + mon_rad / solid_L_vel) + t0 + time_window;

    %Monitoring points for refracted S waves
    mon_end_pts = mon_rad * [sind(refr_S_angle_degs(i)), cosd(refr_S_angle_degs(i))] + [-1;1] * [-cosd(refr_S_angle_degs(i)), sind(refr_S_angle_degs(i))] * src_len / 2;
    tmp = fn_find_nodes_on_line(mod.nds, mon_end_pts(1, :), mon_end_pts(2, :), el_size / 2);
    steps{i}.mon.nds = [steps{i}.mon.nds; tmp; tmp]; %need both DoF at each node in solid
    steps{i}.mon.dfs = [steps{i}.mon.dfs; ones(size(tmp)); ones(size(tmp)) * 2]; %need both DoF at each node in solid
    steps{i}.mon.case = [steps{i}.mon.case; [ones(size(tmp)); ones(size(tmp))] * 4]; %need both DoF at each node in solid
    steps{i}.mon.time_window(4, :) = (src_rad / fluid_vel + mon_rad / solid_L_vel) + t0 + time_window;
end

%Show the mesh
cols = 'cbgm';
if show_mesh
    figure;
    j = 1;
    for i = 1:numel(src_angle_degs)
        display_options.node_sets_to_plot(j).nd = steps{i}.load.frc_nds;
        display_options.node_sets_to_plot(j).col = 'r.';
        j = j + 1;
        display_options.node_sets_to_plot(j).nd = steps{i}.mon.nds(steps{i}.mon.case == 1);
        display_options.node_sets_to_plot(j).col = [cols(j-1), '.'];
        j = j + 1;
        display_options.node_sets_to_plot(j).nd = steps{i}.mon.nds(steps{i}.mon.case == 2);
        display_options.node_sets_to_plot(j).col = [cols(j-1), '.'];
        j = j + 1;
        display_options.node_sets_to_plot(j).nd = steps{i}.mon.nds(steps{i}.mon.case == 3);
        display_options.node_sets_to_plot(j).col = [cols(j-1), '.'];
        j = j + 1;
        display_options.node_sets_to_plot(j).nd = steps{i}.mon.nds(steps{i}.mon.case == 4);
        display_options.node_sets_to_plot(j).col = [cols(j-1), '.'];
        h_patch = fn_show_geometry(mod, matls, display_options);
    end
end

%--------------------------------------------------------------------------
%RUN THE MODEL
res = fn_BristolFE_v2(mod, matls, steps, fe_options);

%--------------------------------------------------------------------------
%SHOW THE RESULTS

%Show the history output as a function of time - here we just sum over all 
%the nodes where displacments were recorded
figure;
for i = 1:numel(src_angle_degs)
    for j = 1:4
        subplot(numel(src_angle_degs), 4, (i - 1) * 4 + j);
        k = steps{i}.mon.case == j;
        d = steps{i}.mon.dfs;
        switch j
            case {1,2}
                tmp = mean(res{i}.dsps(k, :)); %fluid monitoring - no need to worry about DoFs
            case 3
                %solid displacements for L waves - need to resolve onto polarisation direction
                tmp = mean(res{i}.dsps(k & d == 1, :)) * sind(refr_L_angle_degs(i)) + ...
                      mean(res{i}.dsps(k & d == 2, :)) * cosd(refr_L_angle_degs(i));
            case 4
                %solid displacements for S waves - need to resolve onto polarisation direction
                tmp = mean(res{i}.dsps(k & d == 1, :)) * sind(refr_S_angle_degs(i) + 90) + ...
                      mean(res{i}.dsps(k & d == 2, :)) * cosd(refr_S_angle_degs(i) + 90);
        end
        tmp = abs(fn_hilbert(tmp(:)));

        if j == 1
            norm_val = max(abs(tmp));
        end
        tmp = tmp / norm_val;
        val(i, j) = max(tmp(t > steps{i}.mon.time_window(j, 1) & t < steps{i}.mon.time_window(j, 2)));
        semilogy(steps{i}.load.time * 1e6, tmp);
        hold on;
        semilogy([1,1] * steps{i}.mon.time_window(j,1)' * 1e6, [1e-9, 1], 'g')
        semilogy([1,1] * steps{i}.mon.time_window(j,2)' * 1e6, [1e-9, 1], 'g')
        semilogy(steps{i}.mon.time_window(j,:) * 1e6, val(i, j) * [1, 1], 'g');
        axis([0, max(t) * 1e6, 1e-9,1]);
        if i < numel(src_angle_degs)
            set(gca, 'XTickLabel', []);
        else
            xlabel('Time (\mus)')
        end
    end
end

figure;
fine_angle_degs = 0: max(src_angle_degs);
z_fluid = fluid_rho * fluid_vel;
z_solid = solid_rho * solid_L_vel;
% (z_solid - z_fluid) / (z_solid + z_fluid)
col = 'bmg'
for j = 2:4
    semilogy(src_angle_degs, val(:, j), [cols(j), 'o-']);
    hold on;
    ylim([1e-9, 1]);
end
legend('Reflection', 'L refraction', 'S refraction')

%Animate result
if show_mesh
    figure;
    display_options.draw_elements = 0; %makes it easier to see waves if element edges not drawn
    h_patch = fn_show_geometry(mod, matls, display_options);
    anim_options.repeat_n_times = 1;
    anim_options.norm_val = 1000;
    for i = 1:numel(src_angle_degs)
        fn_run_animation(h_patch, res{i}.fld, anim_options);
    end
end