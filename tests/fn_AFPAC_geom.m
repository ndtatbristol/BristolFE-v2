function [model_bdry, water_bdry1, water_bdry2, abs_bdry, cladding_bdry] = fn_AFPAC_geom( ...
    model_length, model_height, ...
    h_wall_thick, v_wall_thick, cladding_thick, abs_layer_thick, ...
    int_rad, r_water_thick)

%Corners of overall domain
model_bdry = [
    0, 0
    0, model_height
    model_length, model_height
    model_length, 0];

%Cladding
cladding_bdry = [
    0, 0
    model_length, 0
    model_length, cladding_thick
    0, cladding_thick];

%Absorbing layer
abs_bdry = [
    abs_layer_thick, 0
    abs_layer_thick, model_bdry(2,2) - abs_layer_thick
    model_bdry(3,1) - abs_layer_thick, model_bdry(2,2) - abs_layer_thick
    model_bdry(3,1) - abs_layer_thick, 0];

%Water (main bit)
water_bdry1 = [
    model_bdry(3,1) - v_wall_thick - int_rad, h_wall_thick
    0, h_wall_thick
    0, model_bdry(3,2)
    model_bdry(3,1) - v_wall_thick, model_bdry(3,2)
    model_bdry(3,1) - v_wall_thick, h_wall_thick + int_rad
    ];
a = linspace(0, -pi /2, 91)';
water_bdry1 = [
    water_bdry1;
    [cos(a), sin(a)] * int_rad + [model_bdry(3,1) - v_wall_thick - int_rad, h_wall_thick + int_rad]
    ];

%Water (RHS)
water_bdry2 = [
    model_bdry(3,1) - r_water_thick, h_wall_thick + int_rad
    model_bdry(3,1) - r_water_thick, model_bdry(3,2)
    model_bdry(3,1), model_bdry(3,2)
    model_bdry(3,1), h_wall_thick
    model_bdry(3,1) - r_water_thick + int_rad, h_wall_thick
    ];
a = linspace(-pi/2, -pi, 91)';
water_bdr2 = [
    water_bdry2;
    [cos(a), sin(a)] * int_rad + [model_bdry(3,1) - r_water_thick + int_rad, h_wall_thick + int_rad]
    ];


end
