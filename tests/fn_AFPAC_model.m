function main = fn_AFPAC_model( ...
    matls, centre_freq, els_per_wavelength, model_length, model_height, ...
    h_wall_thick, v_wall_thick, cladding_thick, abs_layer_thick, ...
    int_rad, r_water_thick, trans_cent, trans_size, trans_angd, subdomain)

do_defect_cases = 1;
do_no_defect_cases = 1;

main.matls = matls;

safety_factor = 3;

%Corners of overall domain
model_corners = [
    0, 0
    0, model_height
    model_length, model_height
    model_length, 0];

%Create root mesh with nodes and materials only
bounding_nds = [min(model_corners); max(model_corners)];



[max_vel, min_vel] = fn_estimate_max_min_vels(matls);

main.mod.design_min_vel = min_vel;
main.mod.design_max_vel = max_vel;
main.mod.design_centre_freq = centre_freq;

%Element size and max time step
main.mod.el_size = min_vel / centre_freq / els_per_wavelength;
main.mod.max_safe_time_step = main.mod.el_size / max_vel / safety_factor;

%Background mesh
[main.mod.nds, main.mod.els] = fn_isometric_structured_mesh(bounding_nds, main.mod.el_size);
main.mod.el_abs_i = zeros(size(main.mod.els, 1), 1);

%First make them all matl 1
main.mod.el_typ_i = repmat({'CPE3'}, [size(main.mod.els, 1), 1]);
main.mod.el_mat_i = ones(size(main.mod.els, 1), 1);

%Then make all the ones below y=cladding_thick matl 2
el_centres = fn_calc_element_centres(main.mod.nds, main.mod.els);
d = el_centres(:,2) < cladding_thick;
main.mod.el_mat_i(d) = 2;

%Absorbing layer
abs_bdry = [
    abs_layer_thick, 0
    abs_layer_thick, model_corners(2,2) - abs_layer_thick
    model_corners(3,1) - abs_layer_thick, model_corners(2,2) - abs_layer_thick
    model_corners(3,1) - abs_layer_thick, 0];
main.mod.el_abs_i = fn_dist_point_to_bdry_2D(el_centres, abs_bdry) / abs_layer_thick;
[d, ~] = fn_elements_in_region(main.mod.nds, main.mod.els, abs_bdry);
main.mod.el_abs_i(d) = 0;
main.mod.el_abs_i(main.mod.el_abs_i > 1) = 1;
main.mod.el_abs_i(main.mod.el_abs_i < 0) = 0;

%Water (main bit)
water_bdry = [
    model_corners(3,1) - v_wall_thick - int_rad, h_wall_thick
    0, h_wall_thick
    0, model_corners(3,2)
    model_corners(3,1) - v_wall_thick, model_corners(3,2)
    model_corners(3,1) - v_wall_thick, h_wall_thick + int_rad
    ];
a = linspace(0, -pi /2, 91)';
water_bdry = [
    water_bdry;
    [cos(a), sin(a)] * int_rad + [model_corners(3,1) - v_wall_thick - int_rad, h_wall_thick + int_rad]
    ];

d = el_centres(:,1) <  model_corners(3,1) - v_wall_thick & el_centres(:,2) > h_wall_thick;
[d, ~] = fn_elements_in_region(main.mod.nds, main.mod.els, water_bdry);
main.mod.el_mat_i(d) = 3;
main.mod.el_typ_i(d) = {'AC2D3'};

%Water (RHS)
water_bdry = [
    model_corners(3,1) - r_water_thick, h_wall_thick + int_rad
    model_corners(3,1) - r_water_thick, model_corners(3,2)
    model_corners(3,1), model_corners(3,2)
    model_corners(3,1), h_wall_thick
    model_corners(3,1) - r_water_thick + int_rad, h_wall_thick
    ];
a = linspace(-pi/2, -pi, 91)';
water_bdry = [
    water_bdry;
    [cos(a), sin(a)] * int_rad + [model_corners(3,1) - r_water_thick + int_rad, h_wall_thick + int_rad]
    ];

d = el_centres(:,1) <  model_corners(3,1) - v_wall_thick & el_centres(:,2) > h_wall_thick;
[d, ~] = fn_elements_in_region(main.mod.nds, main.mod.els, water_bdry);
main.mod.el_mat_i(d) = 3;
main.mod.el_typ_i(d) = {'AC2D3'};

%Next the interface elements
main.mod = fn_add_fluid_solid_interface_els(main.mod);

%Transducer
t = 1;
trans1  = trans_cent - trans_size / 2 * [cosd(trans_angd), sind(trans_angd)];
trans2  = trans_cent + trans_size / 2 * [cosd(trans_angd), sind(trans_angd)];
[main.tx{t}.nds, s] = fn_find_nodes_on_line(main.mod.nds, trans1, trans2, main.mod.el_size / 2);
main.tx{t}.dfs= ones(size(main.tx{t}.nds)) * 4;

% main.steps{t}.load.frc_nds = fn_find_nodes_on_line(main.mod.nds, trans1, trans2, main.mod.el_size / 2);
% main.steps{t}.load.frc_dfs = ones(size(main.steps{t}.load.frc_nds)) * 4;

%Subdomains
a = linspace(0,2*pi, 361)';
for d = 1:numel(subdomain)
    inner_bdry = subdomain(d).cent + [cos(a), sin(a)] * subdomain(d).inner_rad;
    main.doms{d} = fn_create_L_domain(main.mod, inner_bdry, abs_layer_thick);
    % main.doms{d}.mod.bndry_pts = inner_bdry;
    % main.doms{d}.mod.centre = subdomain(d).cent;
end



for d = 1:numel(subdomain)
    s = 0;
    if do_defect_cases
        s = s + 1;
        main.doms{d}.scats{s}.mod.nds = main.doms{d}.mod.nds;
        main.doms{d}.scats{s}.mod.els = main.doms{d}.mod.els;
        main.doms{d}.scats{s}.mod.el_mat_i = main.doms{d}.mod.el_mat_i;
        main.doms{d}.scats{s}.mod.el_abs_i = main.doms{d}.mod.el_abs_i;
        main.doms{d}.scats{s}.mod.el_typ_i = main.doms{d}.mod.el_typ_i;
        main.doms{d}.scats{s}.mod.bdry_lyrs = main.doms{d}.mod.bdry_lyrs;
        main.doms{d}.scats{s}.mod.main_nd_i = main.doms{d}.mod.main_nd_i;
        main.doms{d}.scats{s}.mod.main_el_i = main.doms{d}.mod.main_el_i;

        switch d
            case 1 %surface breaking volumetric defect
                %Remove the existing interface elements
                k = strcmp(main.doms{d}.scats{s}.mod.el_typ_i, 'ASI2D2');
                main.doms{d}.scats{s}.mod.el_typ_i(k) = [];
                main.doms{d}.scats{s}.mod.el_mat_i(k) = [];
                main.doms{d}.scats{s}.mod.el_abs_i(k) = [];
                main.doms{d}.scats{s}.mod.main_el_i(k) = [];
                main.doms{d}.scats{s}.mod.els(k, :) = [];
                %Define defect region
                a = linspace(0, 2*pi, 41)'; a = a(1:end-1);
                r = (0.5 + randn(size(a)) * 0.1) * subdomain(d).inner_rad;
                defect_boundary_points = [r .* cos(a), r .* sin(a)] + subdomain(d).cent;
                el_in_defect = fn_elements_in_region(main.doms{d}.scats{s}.mod.nds, main.doms{d}.scats{s}.mod.els, defect_boundary_points);
                %change material and element type within defect to water
                main.doms{d}.scats{s}.mod.el_mat_i(el_in_defect) = 3;
                main.doms{d}.scats{s}.mod.el_typ_i(el_in_defect) = {'AC2D3'};
                %Add new interface elements
                main.doms{d}.scats{s}.mod = fn_add_fluid_solid_interface_els(main.doms{d}.scats{s}.mod);
            case 2 %Branched crack on cladding side
                % %Define defect region
                % x = [-0.1, -0.2, -0.5, 0.0, 0.4, 0.2, 0.0]' * subdomain(i).inner_rad / sqrt(2) + subdomain(i).cent(1);
                % y = [ 0.0,  0.5,  1.0, 0.6, 0.9, 0.5, 0.0]' * subdomain(i).inner_rad / sqrt(2) + subdomain(i).cent(2);
                % defect_boundary_points = [x, y];
                %Remove elements - THIS NOT WORKING FOR SOME REASON NODES
                %MIXED UP. FILLING WITH AIR INSTEAD
                % el_in_defect = fn_elements_in_region(main.doms{i}.scats{j}.mod.nds, main.doms{i}.scats{j}.mod.els, defect_boundary_points);
                % [~, ~, main.doms{i}.scats{j}.mod.els, main.doms{i}.scats{j}.mod.el_mat_i, main.doms{i}.scats{j}.mod.el_abs_i, main.doms{i}.scats{j}.mod.el_typ_i] = ...
                %     fn_remove_unused_elements(~el_in_defect, main.doms{i}.scats{j}.mod.els, main.doms{i}.scats{j}.mod.el_mat_i, main.doms{i}.scats{j}.mod.el_abs_i, main.doms{i}.scats{j}.mod.el_typ_i);
                % [main.doms{i}.scats{j}.mod.nds, main.doms{i}.scats{j}.mod.els, old_nds, ~] = fn_remove_unused_nodes(main.doms{i}.scats{j}.mod.nds, main.doms{i}.scats{j}.mod.els);
                % main.doms{i}.scats{j}.mod.bdry_lyrs = main.doms{i}.scats{j}.mod.bdry_lyrs(old_nds);
                % main.doms{i}.scats{j}.mod.main_nd_i = main.doms{i}.scats{j}.mod.main_nd_i(old_nds);
                % main.doms{i}.scats{j}.mod.main_el_i = main.doms{i}.scats{j}.mod.main_el_i(old_nds);

                %Remove the existing interface elements
                k = strcmp(main.doms{d}.scats{s}.mod.el_typ_i, 'ASI2D2');
                main.doms{d}.scats{s}.mod.el_typ_i(k) = [];
                main.doms{d}.scats{s}.mod.el_mat_i(k) = [];
                main.doms{d}.scats{s}.mod.el_abs_i(k) = [];
                main.doms{d}.scats{s}.mod.main_el_i(k) = [];
                main.doms{d}.scats{s}.mod.els(k, :) = [];
                %Define defect region
                x = [-0.1, -0.1, -0.5, 0.0, 0.4, 0.1, 0.1]' * subdomain(d).inner_rad / sqrt(2) + subdomain(d).cent(1);
                y = [ 0.0,  0.5,  1.0, 0.6, 0.9, 0.5, 0.0]' * subdomain(d).inner_rad / sqrt(2) + subdomain(d).cent(2);
                defect_boundary_points = [x, y];
                el_in_defect = fn_elements_in_region(main.doms{d}.scats{s}.mod.nds, main.doms{d}.scats{s}.mod.els, defect_boundary_points);
                %change material and element type within defect to water
                main.doms{d}.scats{s}.mod.el_mat_i(el_in_defect) = 4;
                main.doms{d}.scats{s}.mod.el_typ_i(el_in_defect) = {'AC2D3'};
                %Add new interface elements
                % [       main.doms{i}.scats{j}.mod.els, main.doms{i}.scats{j}.mod.el_mat_i, main.doms{i}.scats{j}.mod.el_abs_i, main.doms{i}.scats{j}.mod.el_typ_i] = ...
                %     fn_add_fluid_solid_interface_els(...
                %         main.doms{i}.scats{j}.mod.els, main.doms{i}.scats{j}.mod.el_mat_i, main.doms{i}.scats{j}.mod.el_abs_i, main.doms{i}.scats{j}.mod.el_typ_i, main.doms{i}.scats{j}.mod.nds);
                main.doms{d}.scats{s}.mod = fn_add_fluid_solid_interface_els(main.doms{d}.scats{s}.mod);
            case 3
                %Remove the existing interface elements
                k = strcmp(main.doms{d}.scats{s}.mod.el_typ_i, 'ASI2D2');
                main.doms{d}.scats{s}.mod.el_typ_i(k) = [];
                main.doms{d}.scats{s}.mod.el_mat_i(k) = [];
                main.doms{d}.scats{s}.mod.el_abs_i(k) = [];
                main.doms{d}.scats{s}.mod.main_el_i(k) = [];
                main.doms{d}.scats{s}.mod.els(k, :) = [];
                %Define defect region
                a = linspace(0, 2*pi, 5)'; a = a(1:end-1) + pi / 4;
                r = [0.1, 0.9, 0.1, 0.9]' * subdomain(d).inner_rad;
                defect_boundary_points = [r .* cos(a), r .* sin(a)] + subdomain(d).cent;
                el_in_defect = fn_elements_in_region(main.doms{d}.scats{s}.mod.nds, main.doms{d}.scats{s}.mod.els, defect_boundary_points);
                %change material and element type within defect to water
                main.doms{d}.scats{s}.mod.el_mat_i(el_in_defect) = 3;
                main.doms{d}.scats{s}.mod.el_typ_i(el_in_defect) = {'AC2D3'};
                %Add new interface elements
                % [       main.doms{i}.scats{j}.mod.els, main.doms{i}.scats{j}.mod.el_mat_i, main.doms{i}.scats{j}.mod.el_abs_i, main.doms{i}.scats{j}.mod.el_typ_i] = ...
                %     fn_add_fluid_solid_interface_els(...
                %         main.doms{i}.scats{j}.mod.els, main.doms{i}.scats{j}.mod.el_mat_i, main.doms{i}.scats{j}.mod.el_abs_i, main.doms{i}.scats{j}.mod.el_typ_i, main.doms{i}.scats{j}.mod.nds);
                main.doms{d}.scats{s}.mod = fn_add_fluid_solid_interface_els(main.doms{d}.scats{s}.mod);
        end
        main.doms{d}.scats{s}.mod = fn_tidy_L_model(main.doms{d}.scats{s}.mod);
        % [main.doms{i}.scats{j}.mod.nds, main.doms{i}.scats{j}.mod.els, old_nds, new_nds] = fn_remove_unused_nodes(main.doms{i}.scats{j}.mod.nds, main.doms{i}.scats{j}.mod.els);
        % main.doms{i}.scats{j}.mod.bdry_lyrs = main.doms{i}.mod.bdry_lyrs(old_nds);
        % main.doms{i}.scats{j}.mod.main_nd_i = main.doms{i}.scats{j}.mod.main_nd_i(old_nds);
        % main.doms{i}.scats{j}.mod.main_el_i = main.doms{i}.scats{j}.mod.main_el_i(old_nds);
    end
    if do_no_defect_cases
        s = s + 1;
        %no defects
        main.doms{d}.scats{s}.mod.nds = main.doms{d}.mod.nds;
        main.doms{d}.scats{s}.mod.els = main.doms{d}.mod.els;
        main.doms{d}.scats{s}.mod.el_mat_i = main.doms{d}.mod.el_mat_i;
        main.doms{d}.scats{s}.mod.el_abs_i = main.doms{d}.mod.el_abs_i;
        main.doms{d}.scats{s}.mod.el_typ_i = main.doms{d}.mod.el_typ_i;
        main.doms{d}.scats{s}.mod.bdry_lyrs = main.doms{d}.mod.bdry_lyrs;
        main.doms{d}.scats{s}.mod.main_nd_i = main.doms{d}.mod.main_nd_i;
        main.doms{d}.scats{s}.mod.main_el_i = main.doms{d}.mod.main_el_i;
        % [main.doms{i}.scats{j}.mod.nds, main.doms{i}.scats{j}.mod.els, old_nds, new_nds] = fn_remove_unused_nodes(main.doms{i}.scats{j}.mod.nds, main.doms{i}.scats{j}.mod.els);
        % main.doms{i}.scats{j}.mod.bdry_lyrs = main.doms{i}.mod.bdry_lyrs(old_nds);
        % main.doms{i}.scats{j}.mod.main_nd_i = main.doms{i}.scats{j}.mod.main_nd_i(old_nds);
        % main.doms{i}.scats{j}.mod.main_el_i = main.doms{i}.scats{j}.mod.main_el_i(old_nds);
        main.doms{d}.scats{s}.mod = fn_tidy_L_model(main.doms{d}.scats{s}.mod);
    end
end

end
