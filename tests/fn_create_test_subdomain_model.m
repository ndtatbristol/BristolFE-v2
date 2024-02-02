function main = fn_create_test_subdomain_model(els_per_wavelength, safety_factor, subdoms_to_do)

%Material properties (SI units used throughout)
water_velocity = 1500;
water_density = 1000;
air_velocity = 340;
air_density = 1;

%centre frequency (used for element size calc)
centre_freq = 5e6;
number_of_cycles = 5;


%main geometry - overall bounds
model_height = 20e-3;
model_length = 40e-3;
h_wall_thick = 12e-3;
v_wall_thick = 10e-3;
r_water_thick = 5e-3;
cladding_thick = 1e-3;
int_radius = 3e-3;
abs_bdry_thickness = 2e-3;

%transducer
trans_size = 6e-3;
trans_cent = [abs_bdry_thickness + trans_size / 2, 15e-3];
trans_angd = 20.6; %angle from horisontal: 20.6 for S and 10.8 to L at 45

% time_pts = round(els_per_wavelength / 6 * 20000); %put up to 14k for 6 el/lambda

%Materials
steel_matl_i = 1;
main.matls(steel_matl_i).rho = 8900;
main.matls(steel_matl_i).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3);
main.matls(steel_matl_i).col = hsv2rgb([2/3,0,0.80]);
main.matls(steel_matl_i).name = 'Steel';
main.matls(steel_matl_i).el_typ = 'CPE3';

gold_matl_i = 2;
main.matls(gold_matl_i).rho = 7000;
main.matls(gold_matl_i).D = fn_isotropic_plane_strain_stiffness_matrix(210e9, 0.3);
main.matls(gold_matl_i).col = [0.8672, 0.9375, 0.1875];
main.matls(gold_matl_i).name = 'Gold';
main.matls(gold_matl_i).el_typ = 'CPE3';

water_matl_i = 3;
main.matls(water_matl_i).rho = water_density;
main.matls(water_matl_i).D = water_velocity ^ 2 * water_density;
main.matls(water_matl_i).col = hsv2rgb([0.6,0.5,0.8]);
main.matls(water_matl_i).name = 'Water';
main.matls(water_matl_i).el_typ = 'AC2D3';

% air_matl_i = 4;
% main.matls(air_matl_i).rho = air_density;
% main.matls(air_matl_i).D = air_velocity ^ 2 * water_density;
% main.matls(air_matl_i).col = [1.0, 1.0, 1.0];
% main.matls(air_matl_i).name = 'Air';
% main.matls(air_matl_i).el_typ = 'AC2D3';

%Define sub-domains

%1 At ray entry point
d = 1;
h = trans_cent(2) - h_wall_thick;

if any(strcmp('A', subdoms_to_do))
    % A front wall
    subdomain(d).cent = [trans_cent(1) + h * sind(trans_angd), h_wall_thick];
    subdomain(d).inner_rad = 2e-3;
    d = d + 1;
end

if any(strcmp('B', subdoms_to_do))
    % B back wall
    subdomain(d).cent = [trans_cent(1) + h * sind(trans_angd) + h_wall_thick, 0];
    subdomain(d).inner_rad = 2e-3;
    d = d + 1;
end

if any(strcmp('C', subdoms_to_do))
    % C On radius
    subdomain(d).cent = [trans_cent(1) + h * sind(trans_angd) + 2 * h_wall_thick, h_wall_thick];
    subdomain(d).inner_rad = 2e-3;
    d = d + 1;
end

%--------------------------------------------------------------------------

%Get material bdrys
[model_bdry, water_bdry1, water_bdry2, abs_bdry_pts, cladding_bdry_pts] = fn_AFPAC_geom( ...
    model_length, model_height, ...
    h_wall_thick, v_wall_thick, cladding_thick, abs_bdry_thickness, ...
    int_radius, r_water_thick);

%Work out element size
el_size = fn_get_suitable_el_size(main.matls, centre_freq, els_per_wavelength);

%Create the nodes and elements of the mesh
main.mod = fn_isometric_structured_mesh(model_bdry, el_size);
main.mod.max_safe_time_step = fn_get_suitable_time_step(main.matls, el_size, safety_factor);
main.mod.design_centre_freq = centre_freq;

%First set material of all elements to steel then set elements inside water 
%boundary material to water
main.mod.el_mat_i(:) = steel_matl_i;
main.mod = fn_set_els_inside_bdry_to_mat(main.mod, water_bdry1, water_matl_i);
main.mod = fn_set_els_inside_bdry_to_mat(main.mod, water_bdry2, water_matl_i);
main.mod = fn_set_els_inside_bdry_to_mat(main.mod, cladding_bdry_pts, gold_matl_i);

%Add interface elements - this is crucial otherwise there will be no
%coupling between fluid and solid
main.mod = fn_add_fluid_solid_interface_els(main.mod, main.matls);

%Define the absorbing layer
main.mod = fn_add_absorbing_layer(main.mod, abs_bdry_pts, abs_bdry_thickness);

%Define transducer
no_array_els = 1;
array_el_size = trans_size / no_array_els;
el_cents = linspace(-0.5, 0.5, no_array_els)' * (trans_size - array_el_size) * [cosd(trans_angd), sind(trans_angd)] + trans_cent;
for e = 1:no_array_els
    trans1  = el_cents(e, :) - array_el_size / 2 * [cosd(trans_angd), sind(trans_angd)];
    trans2  = el_cents(e, :) + array_el_size / 2 * [cosd(trans_angd), sind(trans_angd)];
    [main.trans{e}.nds, s] = fn_find_nodes_on_line(main.mod.nds, trans1, trans2, el_size / 2);
    main.trans{e}.dfs = ones(size(main.trans{e}.nds)) * 4;
end

%Add subdomains
a = linspace(0,2*pi,361)';
for d = 1:numel(subdomain)
    inner_bdry = subdomain(d).cent + subdomain(d).inner_rad * [cos(a), sin(a)];
    main.doms{d} = fn_create_subdomain(main.mod, main.matls, inner_bdry, abs_bdry_thickness);
end

end