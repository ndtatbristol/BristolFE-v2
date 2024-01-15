function stiffness_matrix = fn_isotropic_plane_stress_stiffness_matrix(youngs_modulus, poissons_ratio)
%SUMMARY
%   Returns 3 x 3 plane stress stiffness matrix for isotropic material
%   relating stress [sigma_xx, sigma_yy, sigma_xy] to strain [epsilon_xx,
%   epsilon_yy, epsilon_xy].
%INPUTS
%   youngs_modulus - Young's modulus
%   poissons_ratio - Poisson's ratio
%OUTPUTS
%   stiffness_matrix - 3 x 3 stiffness matrix
%--------------------------------------------------------------------------

%input error checks
if ~isscalar(youngs_modulus) | ~isscalar(poissons_ratio)
    error('Youngs modulus and Poissons ratio must be scalars');
end

stiffness_matrix = youngs_modulus / (1 - poissons_ratio ^ 2) * ...
    [1,                     poissons_ratio,         0; ...
    poissons_ratio,         1,                      0; ...
    0,                      0,                      (1 - poissons_ratio) / 2];
end