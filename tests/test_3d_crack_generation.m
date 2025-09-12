clear
close all
addpath(genpath('../code'));

%Simple mesh
model_size = 10;
crnr_pts = [
    0, 0, 0
    1, 1, 1] * model_size;
el_size = 2;
cod = el_size / 10;

steel_mat_i = 1;
matls{steel_mat_i}.rho = 8900; %Density
matls{steel_mat_i}.D = fn_isotropic_stiffness_matrix(210e9, 0.3); 
matls{steel_mat_i}.col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls{steel_mat_i}.name = 'Steel';

el_typ_solid = 'C3D8R'; %C3D8 8 noded brick

mod = fn_3d_structured_mesh_hexahedral_els(crnr_pts, el_size);
mod.el_types = {el_typ_solid};
mod.el_typ_i(:) = find(strcmp(el_typ_solid, mod.el_types));
mod.el_mat_i(:) = steel_mat_i;

%Crack nodes
crack_vtcs1 = [
     1, 2, 3
     8, 9, 6
    -1, 8, 4
    9, 2, 2];
crack_fcs1 = [
    1, 2, 3
    1, 2, 4];

mod = fn_3d_add_crack(mod, crack_vtcs1, crack_fcs1, cod);

figure;
display_options.transparency = 0.5;
display_options.draw_elements = 0;
h_patch = fn_show_geometry(mod, matls, display_options);
patch('Faces', crack_fcs1, 'Vertices', crack_vtcs1,'FaceColor', 'r', 'FaceAlpha', 0.5);

axis on;
xlabel('x')
ylabel('y')
zlabel('z')