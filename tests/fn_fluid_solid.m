function [R_P, T_L, T_S, fluid_theta_ref, solid_theta_L, solid_theta_S] = fn_fluid_solid(fluid_theta_inc, fluid_vel, fluid_rho, solid_vel_L, solid_vel_S, solid_rho)
%Reflection and transmission coefficents for fluid to solid as pressure ratios
%Incident angle in rads between -pi/2 and pi/2
%Schmerr Eq. 6.122
solid_theta_L = asin(sin(fluid_theta_inc) * solid_vel_L / fluid_vel);
solid_theta_S = asin(sin(fluid_theta_inc) * solid_vel_S / fluid_vel);

delta_1 = cos(solid_theta_L);
delta_2 = (solid_vel_L * solid_rho) / (fluid_vel * fluid_rho) * cos(fluid_theta_inc) .* ...
    (cos(2 * solid_theta_S) .^ 2 + solid_vel_S ^ 2 * sin(2 * solid_theta_S) .* sin(2 * solid_theta_L) / solid_vel_L ^ 2);
delta = delta_1 + delta_2;

R_P = (delta_2 - delta_1) ./ delta;

T_L = -2 * (solid_vel_L * solid_rho) * cos(fluid_theta_inc)  .* cos(2 * solid_theta_S) ./ ...
    (fluid_vel * fluid_rho * delta);

T_S = 4 * (solid_vel_S * solid_rho) * cos(fluid_theta_inc) .* cos(solid_theta_L) .* sin(solid_theta_S) ./ ...
    (fluid_vel * fluid_rho * delta);

%all angles are relative to normal pointing into fluid
fluid_theta_ref = -fluid_theta_inc;
solid_theta_L = solid_theta_L + pi;
solid_theta_S = solid_theta_S + pi;
end