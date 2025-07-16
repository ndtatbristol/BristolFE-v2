clear
close all
addpath(genpath('../code'));

%Simple mesh
bdry_pts = [
    0,  0
    0,  10
    10, 10
    10, 0];
el_size = 0.05;
mod = fn_isometric_structured_mesh(bdry_pts, el_size);
steel_mat_i = 1;
matls(steel_mat_i).rho = 8900; %Density
matls(steel_mat_i).D = fn_isotropic_stiffness_matrix(210e9, 0.3); 
matls(steel_mat_i).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(steel_mat_i).name = 'Steel';
matls(steel_mat_i).el_typ = 'C3D8R'; %C3D8 8 noded brick

%Crack nodes
crack_vtcs = [
    2, 3
    4, 4
    7, 8];

[mod, el_cents, ep, en, crack_nds] = fn_add_crack_2d(mod, crack_vtcs, el_size );

crack_vtcs = [
    2, 6
    8, 1];

[mod, el_cents, ep, en, crack_nds] = fn_add_crack_2d(mod, crack_vtcs, el_size );




figure;
options = [];
% options.draw_elements = 1;
fn_show_geometry(mod, matls, options);
hold on;
plot(crack_vtcs(:, 1), crack_vtcs(:, 2), 'r:')
% plot(el_cents(ep, 1), el_cents(ep, 2), 'r+', 'MarkerSize', 3);
% plot(el_cents(en, 1), el_cents(en, 2), 'r_', 'MarkerSize', 3);
% plot(mod.nds(crack_nds, 1), mod.nds(crack_nds, 2), 'ko', 'MarkerSize', 5)
