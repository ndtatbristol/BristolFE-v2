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

%Defined crack vertices - a nice branched crack with some roughness
pts1 = 20;
len1 = 5;
pts2 = 10;
len2 = 2.5;
%top part
crack_vtcs1 = [-el_size, 5] + fn_2d_random_walk(pts1, len1/pts1, 0, linspace(0,1,pts1) * pi/8, 0.4);
%branch start (mid way along top branch)
crack_vtcs2 = crack_vtcs1(round(pts1 / 2), :) + fn_2d_random_walk(pts2, len2/pts2, 0, -pi/8-linspace(0,1,pts2) * pi/8, 0.4);

%Create the mesh
mod = fn_2d_structured_mesh_triangular_els(bdry_pts, el_size);
mod.el_types = {el_typ_solid};
mod.el_typ_i(:) = find(strcmp(el_typ_solid, mod.el_types));
mod.el_mat_i(:) = steel_mat_i;

%Add the cracks
mod = fn_2d_add_crack(mod, crack_vtcs1, [], cod);
mod = fn_2d_add_crack(mod, crack_vtcs2, [], cod);

%Plot result
figure;
options = [];
fn_show_geometry(mod, matls, options);
hold on;
plot(crack_vtcs1(:, 1), crack_vtcs1(:, 2), 'r:')
plot(crack_vtcs2(:, 1), crack_vtcs2(:, 2), 'r:')
