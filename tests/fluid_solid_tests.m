%Fluid-solid transmission, reflection test
clear;
close all;
restoredefaultpath;
addpath(genpath('../code'));

%--------------------------------------------------------------------------
%GENERAL PARAMETERS

%Incident angles (all angles are measured relative to normal pointing into
%fluid)
src_angle_degs = [0:5:15]; %0 is normal incidence
% src_angle_degs = [0: 0.5: 30];
src_angle_degs = 14;

abs_bdry_thickness = 1e-3;
model_size = 5e-3;

%Source transducer radial position
src_rad = model_size - 2 * abs_bdry_thickness;

%Source transducer length
src_trans_len = model_size * 1;

%Monitoring transducer radial position
mon_rad = model_size - 3 * abs_bdry_thickness;
mon_trans_len = model_size * 0.5;

centre_freq = 5e6;
no_cycles = 4;
els_per_wavelength = 10;
fe_options.field_output_every_n_frames = 20;
cols = 'cbgmk';

%Material properties
fluid_vel = 1500;
fluid_rho = 1000;
solid_vel_L = 6000;
solid_vel_S = 3000;
solid_rho = 9000;

fluid_theta_inc = src_angle_degs * pi  / 180;
[~, ~, ~, fluid_theta_ref, solid_theta_L, solid_theta_S] = fn_fluid_solid(fluid_theta_inc, fluid_vel, fluid_rho, solid_vel_L, solid_vel_S, solid_rho);

show_mesh_only = 0;
show_animation = 1;

fe_options.interface_damping_factor = 0;
fe_options.field_output_type = 'mean(u1)';
safety_factor = 1.5;


%--------------------------------------------------------------------------
%Prepare materials

steel_matl_i = 1;
matls(steel_matl_i).rho = 8900; %Density
[lambda, mu] = fn_lame_from_velocities_and_density(solid_vel_L, solid_vel_S, solid_rho);
[youngs_modulus, poissons_ratio] = fn_youngs_from_lame(lambda, mu);
matls(steel_matl_i).D = fn_isotropic_stiffness_matrix(youngs_modulus, poissons_ratio); 
matls(steel_matl_i).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(steel_matl_i).name = 'Solid';
matls(steel_matl_i).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

water_matl_i = 2;
matls(water_matl_i).rho = fluid_rho;
matls(water_matl_i).D = fluid_vel ^ 2 * 1000;
matls(water_matl_i).col = hsv2rgb([0.6,0.5,0.8]);
matls(water_matl_i).name = 'Fluid'; 
matls(water_matl_i).el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid

%--------------------------------------------------------------------------

%PREPARE THE MESH

%absorbing boundary
a = linspace(0, 2 * pi, 361)';
abs_bdry_pts = [cos(a), sin(a)] * (model_size - abs_bdry_thickness);

%Work out element size
el_size = fn_get_suitable_el_size(matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh
%Define shape of model
bdry_pts = [cos(a), sin(a)] * model_size;
mod = fn_isometric_structured_mesh(bdry_pts, el_size);

%First set material of all elements to steel then set elements inside water 
%boundary material to water
mod.el_mat_i(:) = steel_matl_i;
water_bdry_pts = [-1, 0; 1, 0; 1, -1; -1, -1] * model_size;
mod = fn_set_els_inside_bdry_to_mat(mod, water_bdry_pts, water_matl_i);

%Add interface elements - this is crucial otherwise there will be no
%coupling between fluid and solid
mod = fn_add_fluid_solid_interface_els(mod, matls);

%Define the absorbing layer
mod = fn_add_absorbing_layer(mod, abs_bdry_pts, abs_bdry_thickness);

%Identify nodes along the source line to say where the loading will be 
%when FE model is run
time_step = fn_get_suitable_time_step(matls, el_size, safety_factor);
max_time = model_size * 2 / fluid_vel;
t = 0: time_step:  max_time;
in_sig = fn_gaussian_pulse(t, centre_freq, no_cycles);
[~, i] = max(abs(fn_hilbert(in_sig(:))));
t0 = t(i);
time_window = [-1, 1] * no_cycles / centre_freq; %time window used for looking for peak amplitudes of reflected/refracted waves


%Some inline functions for calculating useful coordinates
fn_end_pts = @(rad, ang, len)   [sin(ang), -cos(ang)] * rad + [-1; 1] * [cos(ang), sin(ang)] * len / 2;
fn_moni_vec_L = @(ang)          [sin(ang), -cos(ang)]; 
fn_moni_vec_S = @(ang)          [cos(ang),  sin(ang)]; 

src_dir = 4; %direction of forces applied: 1 = x, 2 = y, 3 = z (for solids), 4 = volumetric expansion (for fluids)

for i = 1:numel(src_angle_degs)
    % src_end_pts = src_rad * [sind(src_angle_degs(i)), cosd(src_angle_degs(i))] + [-1;1] * [-cosd(src_angle_degs(i)), sind(src_angle_degs(i))] * src_len / 2;
    src_end_pts = fn_end_pts(src_rad, fluid_theta_inc(i), src_trans_len);
    steps{i}.load.frc_nds = fn_find_nodes_on_line(mod.nds, src_end_pts(1, :), src_end_pts(2, :), el_size / 2);
    steps{i}.load.frc_dfs = ones(size(steps{i}.load.frc_nds)) * src_dir;
    steps{i}.load.time = t;
    steps{i}.load.frcs = in_sig;

    %Now add the monitoring nodes for each sort of wave
    steps{i}.mon.nds = [];
    steps{i}.mon.dfs = [];
    steps{i}.mon.case = [];

    %Monitoring points for incident wave in water
    mon_end_pts = fn_end_pts(mon_rad, fluid_theta_inc(i), mon_trans_len);
    tmp = fn_find_nodes_on_line(mod.nds, mon_end_pts(1, :), mon_end_pts(2, :), el_size / 2);
    steps{i}.mon.nds = [steps{i}.mon.nds; tmp];
    steps{i}.mon.dfs = [steps{i}.mon.dfs; ones(size(tmp)) * src_dir];
    steps{i}.mon.case = [steps{i}.mon.case; ones(size(tmp)) * 1];
    steps{i}.mon.time_window(1, :) = (src_rad - mon_rad) / fluid_vel + t0 + time_window;

    %Monitoring points for reflected waves
    mon_end_pts = fn_end_pts(mon_rad, fluid_theta_ref(i), mon_trans_len);
    tmp = fn_find_nodes_on_line(mod.nds, mon_end_pts(1, :), mon_end_pts(2, :), el_size / 2);
    steps{i}.mon.nds = [steps{i}.mon.nds; tmp];
    steps{i}.mon.dfs = [steps{i}.mon.dfs; ones(size(tmp)) * src_dir];
    steps{i}.mon.case = [steps{i}.mon.case; ones(size(tmp)) * 2];
    steps{i}.mon.time_window(2, :) = (src_rad + mon_rad) / fluid_vel + t0 + time_window;

    %Monitoring points for refracted L waves
    mon_end_pts = fn_end_pts(mon_rad, solid_theta_L(i), mon_trans_len);
    tmp = fn_find_nodes_on_line(mod.nds, mon_end_pts(1, :), mon_end_pts(2, :), el_size / 2);
    steps{i}.mon.nds = [steps{i}.mon.nds; tmp; tmp]; %need both DoF at each node in solid
    steps{i}.mon.dfs = [steps{i}.mon.dfs; ones(size(tmp)); ones(size(tmp)) * 2]; %need both DoF at each node in solid
    steps{i}.mon.case = [steps{i}.mon.case; [ones(size(tmp)); ones(size(tmp))] * 3]; %need both DoF at each node in solid
    steps{i}.mon.time_window(3, :) = (src_rad / fluid_vel + mon_rad / solid_vel_L) + t0 + time_window;

    %Monitoring points for refracted S waves
    mon_end_pts = fn_end_pts(mon_rad, solid_theta_S(i), mon_trans_len);
    tmp = fn_find_nodes_on_line(mod.nds, mon_end_pts(1, :), mon_end_pts(2, :), el_size / 2);
    steps{i}.mon.nds = [steps{i}.mon.nds; tmp; tmp]; %need both DoF at each node in solid
    steps{i}.mon.dfs = [steps{i}.mon.dfs; ones(size(tmp)); ones(size(tmp)) * 2]; %need both DoF at each node in solid
    steps{i}.mon.case = [steps{i}.mon.case; [ones(size(tmp)); ones(size(tmp))] * 4]; %need both DoF at each node in solid
    steps{i}.mon.time_window(4, :) = (src_rad / fluid_vel + mon_rad / solid_vel_S) + t0 + time_window;
end

%Node sets to plot
display_options = [];
j = 1;
for i = 1:numel(src_angle_degs)
    display_options.node_sets_to_plot(j).nd = steps{i}.load.frc_nds;
    display_options.node_sets_to_plot(j).col = 'r.';
    j = j + 1;
    display_options.node_sets_to_plot(j).nd = steps{i}.mon.nds(steps{i}.mon.case == 1);
    display_options.node_sets_to_plot(j).col = [cols(1), '.'];
    j = j + 1;
    display_options.node_sets_to_plot(j).nd = steps{i}.mon.nds(steps{i}.mon.case == 2);
    display_options.node_sets_to_plot(j).col = [cols(2), '.'];
    j = j + 1;
    display_options.node_sets_to_plot(j).nd = steps{i}.mon.nds(steps{i}.mon.case == 3);
    display_options.node_sets_to_plot(j).col = [cols(3), '.'];
    j = j + 1;
    display_options.node_sets_to_plot(j).nd = steps{i}.mon.nds(steps{i}.mon.case == 4);
    display_options.node_sets_to_plot(j).col = [cols(4), '.'];
    j = j + 1;
end

%Show the mesh
if show_mesh_only
    figure;
    h_patch = fn_show_geometry(mod, matls, display_options);
    return
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
        cs = steps{i}.mon.case(res{i}.valid_mon_dsps); %cases associated with actual displacement outputs (nesc in case not all requested node/dof combos are valid
        dof = steps{i}.mon.dfs(res{i}.valid_mon_dsps);
        switch j
            case {1,2}
                tmp = mean(res{i}.dsps(cs == j, :)); %fluid monitoring - no need to worry about DoFs
            case 3
                %solid displacements for L waves - need to resolve onto polarisation direction
                moni_vec_L = fn_moni_vec_L(solid_theta_L(i));
                tmp = res{i}.dsps(cs == j & dof == 1, :) * moni_vec_L(1) + ...
                      res{i}.dsps(cs == j & dof == 2, :) * moni_vec_L(2);
                tmp = mean(tmp) * solid_rho * solid_vel_L;
            case 4
                %solid displacements for S waves - need to resolve onto polarisation direction
                moni_vec_S = fn_moni_vec_S(solid_theta_S(i));
                tmp = res{i}.dsps(cs == j & dof == 1, :) * moni_vec_S(1) + ...
                      res{i}.dsps(cs == j & dof == 2, :) * moni_vec_S(2);
                tmp = mean(tmp) * solid_rho * solid_vel_S;
        end
        tmp = abs(fn_hilbert(tmp(:)));
        %At this point, tmp has units of pressure or stress integrated
        %w.r.t. time in both solid and fluid, so amplitude ratios should be
        %equivalent to pressure transmission and reflection coefficients

        if j == 1
            norm_val = max(abs(tmp));
        end
        tmp = tmp / norm_val;
        fe_coeffs_pressure(i, j) = max(tmp(t > steps{i}.mon.time_window(j, 1) & t < steps{i}.mon.time_window(j, 2)));
        plot(steps{i}.load.time * 1e6, tmp)
        hold on;
        plot([1,1] * steps{i}.mon.time_window(j,1)' * 1e6, [0, 2], 'g')
        plot([1,1] * steps{i}.mon.time_window(j,2)' * 1e6, [0, 2], 'g')
        plot(steps{i}.mon.time_window(j,:) * 1e6, fe_coeffs_pressure(i, j) * [1, 1], 'g');
        if i < numel(src_angle_degs)
            set(gca, 'XTickLabel', []);
        else
            xlabel('Time (\mus)')
        end
    end
end

%Pressure coefficients
figure;
fluid_theta_inc_fine = linspace(fluid_theta_inc(1), fluid_theta_inc(end), 200);
[R_P, T_L, T_S] = fn_fluid_solid(fluid_theta_inc_fine, fluid_vel, fluid_rho, solid_vel_L, solid_vel_S, solid_rho);
theory_coeffs_pressure = [ones(size(R_P)); R_P; T_L; T_S] .';
for j = 2:4
    plot(src_angle_degs, fe_coeffs_pressure(:, j), [cols(j), 'x']);
    hold on;
    plot(fluid_theta_inc_fine * 180 / pi, abs(theory_coeffs_pressure(:, j)), cols(j));
end
xlabel('Angle (^o)')
ylabel('Relf. or refr. coefficient (pressure)')
legend('Reflection', '', 'L refraction', '', 'S refraction', '')

%Power coefficients
figure;
z = [fluid_vel * fluid_rho, fluid_vel * fluid_rho, solid_vel_L * solid_rho, solid_vel_S * solid_rho];
fe_coeffs_power = fe_coeffs_pressure .^ 2 ./ z / 2;
fe_coeffs_power(:, 2:end) = fe_coeffs_power(:, 2:end) ./ fe_coeffs_power(:, 1);
fe_coeffs_power = [fe_coeffs_power, sum(fe_coeffs_power, 2)];
theory_coeffs_power = theory_coeffs_pressure .^ 2 ./ z / 2;
theory_coeffs_power(:, 2:end) = theory_coeffs_power(:, 2:end) ./ theory_coeffs_power(:, 1);
theory_coeffs_power = [theory_coeffs_power, sum(theory_coeffs_power, 2)];
for j = 2:5
    plot(src_angle_degs, fe_coeffs_power(:, j), [cols(j), 'x']);
    hold on;
    plot(fluid_theta_inc_fine * 180 / pi, abs(theory_coeffs_power(:, j)), cols(j));
end
% ylim([0, 1]);
xlabel('Angle (^o)')
ylabel('Relf. or refr. coefficient (power)')
legend('Reflection', '', 'L refraction', '', 'S refraction', '')


%Animate result
if show_animation
    figure;
    display_options.draw_elements = 0; %makes it easier to see waves if element edges not drawn
    h_patch = fn_show_geometry(mod, matls, display_options);
    anim_options = [];
    anim_options.repeat_n_times = 1;
    % anim_options.db_or_linear = 'linear';
    for i = 1:numel(src_angle_degs)
        fn_run_animation(h_patch, res{i}.fld, anim_options);
    end
end

% save('fluid-solid result 20els per wavelength','fe_coeffs_pressure')