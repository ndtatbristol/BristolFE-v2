clear
close all

el_size = 1e-3;

%Following approx for steel water
sol_rho = 8.9e3;
sol_mod = sol_rho * 6e3 ^ 2;
sol_nu = 0.3; %used for the 2D elements in third example
sol_damping = 0;

liq_rho = 1e3;
liq_mod = liq_rho * 1.5e3 ^ 2;
liq_damping = 0;

coupling_eff = 1;

no_time_steps = 101;
min_decades = -3;
max_decades = 1;
ymax = 1.01;

rayleigh_damping = 1000000;

%--------------------------------------------------------------------------

%Critical timestep calc (usual method - min time to transit element) - used
%to determine range of timesteps to investigate
w1 = sqrt(sol_mod / sol_rho);
w2 = sqrt(liq_mod / liq_rho);
crit_dt = min(el_size ./ [w1, w2]) * sqrt(2);
dt = logspace(min_decades,max_decades,no_time_steps)' * crit_dt;

%--------------------------------------------------------------------------
%First test - two-element structure with pinned end-nodes, so only 2DoF are
%at centre node where time-integrated pressure and displacement are coupled
%so system matrices are 2x2.

%Values for element stiffness and mass matrices
k1 = sol_mod / el_size;
m1 = sol_rho * el_size / 2;
c1 = sol_damping;

med2fac = -1/liq_rho;
k2 = liq_mod / el_size * med2fac;
m2 = liq_rho * el_size / 2 * med2fac;
c2 = liq_damping * med2fac;

%Associated global matrices
M = [m1, 0;
     0,  m2];
K = [k1, 0;
     0,  k2];
C = [c1, coupling_eff;
     coupling_eff,  c2];
% uncoupled_C =  [  c1, 0;
%                   0,  c2];

%Uncoupled case first to show critical time step
uncoupled = zeros(numel(dt), 1);
coupled = zeros(numel(dt), 1);

for i = 1:numel(dt)
    [~, uncoupled(i)] = fn_amplification_matrix(M, C .* eye(size(C)), K, dt(i), []);
    [~, coupled(i)] = fn_amplification_matrix(M, C, K, dt(i), []);
end

% figure
% loglog(dt / crit_dt, uncoupled, 'b'); hold on;
% loglog(dt / crit_dt, coupled, 'r.'); hold on;
% ylim([1, ymax]);
% legend({'Steel', 'Coupled'});
% xlabel('\delta / \delta_{crit}')
% ylabel('Max Eigenvalue')
% title('Two 1D elements, pinned ends')

%--------------------------------------------------------------------------
%Second test - as before but with free end-nodes, so now 4DoF totol

%Associated global matrices
M = [m1, 0, 0, 0
     0,  m1, 0, 0
     0,  0, m2, 0
     0,  0, 0, m2];
K = [k1, 0, 0, 0
     0,  k1, 0, 0
     0,  0, k2, 0
     0,  0, 0, k2];
C = [c1, 0, 0, 0
     0,  c1, coupling_eff, 0
     0,  coupling_eff, c2, 0
     0,  0, 0, c2];
% uncoupled_C= [c1, 0, 0, 0
%      0,  c1, 0, 0
%      0,  0, c2, 0
%      0,  0, 0, c2];

%Proposed interface
uncoupled2 = zeros(size(dt));
coupled2 = zeros(size(dt));
for i = 1:numel(dt)
    [~, uncoupled2(i)] = fn_amplification_matrix(M, C .* eye(size(C)), K, dt(i), []);
    [~, coupled2(i)] = fn_amplification_matrix(M, C, K, dt(i), []);
end

% figure
% loglog(dt / crit_dt, uncoupled2, 'b'); hold on;
% loglog(dt / crit_dt, coupled2, 'r.'); hold on;
% ylim([1, ymax]);
% legend({'Steel', 'Coupled'});
% xlabel('\delta / \delta_{crit}')
% ylabel('Max Eigenvalue')
% title('Two 1D elements, free ends')

%--------------------------------------------------------------------------
%Third test - two 2D elements

matls(1).rho = sol_rho; %Density
matls(1).D = fn_isotropic_stiffness_matrix(sol_mod, sol_nu); 
matls(1).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(1).name = 'Steel';
matls(1).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

matls(2).rho = liq_rho;
matls(2).D = liq_mod;
matls(2).col = hsv2rgb([0.6,0.5,0.8]);
matls(2).name = 'Water';
matls(2).el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid

%pair of elements with interface
mod.nds = [0.5, sqrt(3) / 2; 0, 0; 1, 0; 0.5, -sqrt(3)/2] * el_size * 2 * sqrt(3/2);
mod.els = [1,2,3; 2,4,3];
mod.el_mat_i = [1;2];
mod = fn_add_fluid_solid_interface_els(mod, matls, []);

maxEV = zeros(numel(dt), 2);

%Stability vs. time step
for i = 1:numel(dt)
    [K, C, M, gl_lookup] = fn_build_global_matrices_v5(mod.nds, mod.els, mod.el_mat_i, zeros(size(mod.el_mat_i)), mod.el_typ_i, matls, []);
    [~, maxEV(i, 1)] = fn_amplification_matrix(M, C .* eye(size(C)), K, dt(i), []);
    % if rayleigh_damping
    %     C = C + M * rayleigh_damping;
    % end
    C = fn_add_interface_damping(mod, M, C, gl_lookup, rayleigh_damping);
    [~, maxEV(i, 2)] = fn_amplification_matrix(M, C, K, dt(i), []);
end

figure;
subplot(1,2,1)
h_patch = fn_show_geometry(mod, matls, []);
subplot(1,2,2)
loglog(dt / crit_dt, maxEV(:, 1), 'b'); hold on;
loglog(dt / crit_dt, maxEV(:, 2), 'r.');
ylim([1, ymax]);
xlabel('dt / dt_{crit}')
ylabel('Max EV');
legend('Steel', 'Coupled');
title('Two 2D elements')

%--------------------------------------------------------------------------
%Fourth test - a bigger model
clear('mod');
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

mod = fn_isometric_structured_mesh(bdry_pts, el_size * 2 * sqrt(3/2));

%First set material of all elements to steel then set elements inside water 
%boundary material to water
mod.el_mat_i(:) = 1;
mod = fn_set_els_inside_bdry_to_mat(mod, water_bdry_pts, 2);

%Add interface elements - this is crucial otherwise there will be no
%coupling between fluid and solid
mod = fn_add_fluid_solid_interface_els(mod, matls);

maxEV = zeros(numel(dt), 2);
%Stability vs. time step
for i = 1:numel(dt)
    [K, C, M, gl_lookup] = fn_build_global_matrices_v5(mod.nds, mod.els, mod.el_mat_i, zeros(size(mod.el_mat_i)), mod.el_typ_i, matls, []);
    [~, maxEV(i, 1)] = fn_amplification_matrix(M, C .* eye(size(C)), K, dt(i), []);
    % if rayleigh_damping
    %     C = C + M * rayleigh_damping;
    % end
    C = fn_add_interface_damping(mod, M, C, gl_lookup, rayleigh_damping);
    [~, maxEV(i, 2)] = fn_amplification_matrix(M, C, K, dt(i), []);
end

figure;
subplot(1,2,1)
h_patch = fn_show_geometry(mod, matls, []);
subplot(1,2,2)
loglog(dt / crit_dt, maxEV(:, 1), 'b'); hold on;
loglog(dt / crit_dt, maxEV(:, 2), 'r.');
ylim([1, ymax]);
xlabel('dt / dt_{crit}')
ylabel('Max EV');
legend('Steel', 'Coupled');
title('Proper mesh')

