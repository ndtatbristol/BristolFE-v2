function D = fn_isotropic_stiffness_matrix(E, v)
%SUMMARY
%   Returns 6 x 6 stiffness matrix for isotropic material
%   relating stress to strain.
%INPUTS
%   E - Young's modulus
%   v - Poisson's ratio
%OUTPUTS
%   stiffness_matrix - 6 x 6 stiffness matrix
%--------------------------------------------------------------------------

%input error checks
if ~isscalar(E) | ~isscalar(v)
    error('Youngs modulus and Poissons ratio must be scalars');
end

D = E / (1 - 2 * v) / (1 + v) * ...
    [1 - v,     v,          v,          0,               0,               0
    v,          1 - v,      v,          0,               0,               0
    v,          v,          1 - v,      0,               0,               0
    0,          0,          0,          (1 - 2 * v) / 2, 0,               0
    0,          0,          0,          0,               (1 - 2 * v) / 2, 0
    0,          0,          0,          0,               0,               (1 - 2 * v) / 2];

end