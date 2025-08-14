clear
close all
addpath(genpath('../code'));

%Simple mesh
model_size = 10;
crnr_pts = [
    0, 0, 0
    1, 1, 1] * model_size;
el_size = 0.5;
cod = el_size / 10;
mod = fn_3d_cubic_structured_mesh(crnr_pts, el_size);

steel_mat_i = 1;
matls(steel_mat_i).rho = 8900; %Density
matls(steel_mat_i).D = fn_isotropic_stiffness_matrix(210e9, 0.3); 
matls(steel_mat_i).col = hsv2rgb([2/3,0,0.80]); %Colour for display
matls(steel_mat_i).name = 'Steel';
matls(steel_mat_i).el_typ = 'C3D8R'; %C3D8 8 noded brick


figure;
display_options.el_typ_i = {matls(mod.el_mat_i).el_typ}'; %this is irritating - there is no el_typ_i in mod - el_types only accessed through material
display_options.transparency = 0.5;
display_options.draw_elements = 1;
h_patch = fn_show_geometry(mod, matls, display_options);




%Crack nodes
crack_vtcs1 = [
    1, 2, 3
    8, 9, 6
    1, 8, 4];
crack_fcs1 = [
    1, 2, 3];

patch('Faces', crack_fcs1, 'Vertices', crack_vtcs1,'FaceColor', 'r', 'FaceAlpha', 0.5);

[mod, el_cents, ep, en, crack_nds] = fn_add_crack_3d(mod, crack_vtcs1, crack_fcs1, cod);
hold on; plot3(el_cents(ep, 1), el_cents(ep, 2), el_cents(ep, 3), 'r.');
hold on; plot3(el_cents(en, 1), el_cents(en, 2), el_cents(en, 3), 'g.');

return

[mod, el_cents, ep, en, crack_nds] = fn_add_crack_2d(mod, crack_vtcs1, cod);

crack_vtcs2 = [
    2, 6
    8, 1];

[mod, el_cents, ep, en, crack_nds] = fn_add_crack_2d(mod, crack_vtcs2, cod);




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
