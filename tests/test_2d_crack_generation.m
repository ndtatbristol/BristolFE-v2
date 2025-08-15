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
cod = el_size / 10;


steel_mat_i = 1;
matls{steel_mat_i}.rho = 8900; %Density
matls{steel_mat_i}.D = fn_isotropic_stiffness_matrix(210e9, 0.3); 
matls{steel_mat_i}.col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls{steel_mat_i}.name = 'Steel';

el_typ_solid = 'CPE3'; 


%Crack nodes
crack_vtcs1 = [
    2, 3
    4, 4
    7, 8];
crack_vtcs2 = [
    2, 6
    8, 1];
crack_vtcs3 = [
    1, 1
    1, 2
    2,2
    2,1
    1,1];


mod = fn_2d_isometric_structured_mesh(bdry_pts, el_size);
mod.el_types = {el_typ_solid};
mod.el_typ_i(:) = find(strcmp(el_typ_solid, mod.el_types));
mod.el_mat_i(:) = steel_mat_i;

[mod, el_cents, ep, en, crack_nds] = fn_add_crack(mod, crack_vtcs1, [], cod);
[mod, el_cents, ep, en, crack_nds] = fn_add_crack(mod, crack_vtcs2, [], cod);
[mod, el_cents, ep, en, crack_nds] = fn_add_crack(mod, crack_vtcs3, [], cod);




figure;
options = [];
% options.draw_elements = 1;
fn_show_geometry(mod, matls, options);
hold on;
plot(crack_vtcs1(:, 1), crack_vtcs1(:, 2), 'r:')
plot(crack_vtcs2(:, 1), crack_vtcs2(:, 2), 'r:')
% plot(el_cents(ep, 1), el_cents(ep, 2), 'r+', 'MarkerSize', 3);
% plot(el_cents(en, 1), el_cents(en, 2), 'r_', 'MarkerSize', 3);
% plot(mod.nds(crack_nds, 1), mod.nds(crack_nds, 2), 'ko', 'MarkerSize', 5)
