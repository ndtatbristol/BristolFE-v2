clear
% close all
restoredefaultpath;
addpath(genpath('../code'));

%Solid props
sol_rho = 8900;
sol_E = 210e9;
sol_nu = 0.3;

%Liquid props
liq_rho = 1000;
liq_B = 1500^2 * 1000;

el_size = 1e-3;

fe_options.global_matrix_builder_version = 'v5';
fe_options.dynamic_solver_version = 'v6';

fe_options.solver_mode = 'vel at last half time step';
% fe_options.solver_mode = 'vel at curent time step';

%--------------------------------------------------------------------------

%Stability investigation for interface elements
matls(1).rho = sol_rho; %Density
matls(1).D = fn_isotropic_stiffness_matrix(sol_E, sol_nu); 
matls(1).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(1).name = 'Steel';
matls(1).el_typ = 'CPE3'; %CPE3 must be the element type for a solid

matls(2).rho = liq_rho;
matls(2).D = liq_B;
matls(2).col = hsv2rgb([0.6,0.5,0.8]);
matls(2).name = 'Water';
matls(2).el_typ = 'AC2D3'; %AC2D3 must be the element type for a fluid

vmax = sqrt(sol_E .* (1 - sol_nu) ./ sol_rho ./ (1 + sol_nu) ./ (1 - 2 * sol_nu));

%single element of steel
mod{1}.nds = [0.5, sqrt(3) / 2; 0, 0; 1, 0] * el_size;
mod{1}.els = [1,2,3];
mod{1}.el_mat_i = [1];
mod{1}.el_typ_i = {matls(mod{1}.el_mat_i).el_typ};

%single element of water
mod{2}.nds = [0, 0; 1, 0; 0.5, -sqrt(3)/2] * el_size;
mod{2}.els = [1, 2, 3];
mod{2}.el_mat_i = [2];
mod{2}.el_typ_i = {matls(mod{2}.el_mat_i).el_typ};


%pair of elements with interface
mod{3}.nds = [0.5, sqrt(3) / 2; 0, 0; 1, 0; 0.5, -sqrt(3)/2] * el_size;
mod{3}.els = [1,2,3; 2,3,4];
mod{3}.el_mat_i = [1;2];
options.fluid_el_names = {'AC2D3', 'AC2D3_v2'};
mod{3} = fn_add_fluid_solid_interface_els(mod{3}, matls, options);


crit_dt = el_size / vmax;
dt = logspace(-1,1,101)' * crit_dt;
maxEV = zeros(numel(dt), numel(mod));

%Stability vs. time step
for i = 1:numel(dt)
    for j = 1:numel(mod)
        [K, C, M, ~] = fn_build_global_matrices_v5(mod{j}.nds, mod{j}.els, mod{j}.el_mat_i, zeros(size(mod{j}.el_mat_i)), mod{j}.el_typ_i, matls, fe_options);
        [~, maxEV(i, j)] = fn_amplification_matrix(M, C, K, dt(i), fe_options);
    end
end

figure;
loglog(dt / crit_dt, abs(maxEV));
xlabel('dt / dt_{crit}')
ylabel('Max EV');
legend('Steel', 'Water', 'Combined');




% %Stability vs. safety factor
% ev = zeros(numel(safety_factors), numel(mod));
% for i = 1:numel(safety_factors)
%     for j = 1:numel(mod)
%         [K, C, M, ~] = fn_build_global_matrices_v5(mod{j}.nds, mod{j}.els, mod{j}.el_mat_i, zeros(size(mod{j}.el_mat_i)), mod{j}.el_typ_i, matls, fe_options);
%         [~, maxEV(i, j)] = fn_amplification_matrix(M, C, K, dt * safety_factors(i), fe_options);
%     end
% end
% 
% figure;
% for j = 1:numel(mod)
%     subplot(numel(mod), 1, j);
%     semilogy(safety_factors(:), abs(maxEV(:, j)));
%     xlabel('Safety factor');
%     ylabel('Max EV')
% end


return

%Actual runs using params
steps{1}.load.time = [0:1000] * crit_dt;
steps{1}.load.frc_nds = 1;
steps{1}.load.frc_dfs = 1;
steps{1}.load.frcs = zeros(size(steps{1}.load.time));
steps{1}.load.frcs(10) = 1;
steps{1}.mon.nds = 1;
steps{1}.mon.dfs = 1;

for j = 1:numel(mod)
    [res{j}, mats] = fn_BristolFE_v2(mod{j}, matls, steps, fe_options);
    [~, maxEV_sf1(j)] = fn_amplification_matrix(mats.M, mats.C, mats.K, crit_dt, fe_options);
end

figure;
for j = 1:numel(mod)
    subplot(1, numel(mod), j);
    fn_show_geometry(mod{j}, matls, []);
end

figure;
for j = 1:numel(mod)
    subplot(numel(mod), 1, j);
    semilogy(steps{1}.load.time, abs(res{j}{1}.dsps));
    title(sprintf('Max EV = %.8f', maxEV_sf1(j)))
end


