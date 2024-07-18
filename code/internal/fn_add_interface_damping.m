function C = fn_add_interface_damping(els, el_typ_i, K, C, M, gl_lookup, varargin)
%SUMMARY
%   Adds Rayleigh damping to interface nodes in model
%USAGE
%   mats = fn_add_interface_damping(mod, mats [, factor])
%INPUTS
%   els - m x max_nds_per_el matrix of element nodes. The row number is the element
%   number; columns are the node numbers of the nodes for each element with
%   trailing zeros if element has less nodes than max_nds_per_el
%   el_typ_i - m x 1 cell array of element types
%   K, C, M - global stiffness, mass and damping matrices
%   gl_lookup - n_nds x n_dof matrix of global matrix indices for node-DoF
%   pairs
%   [factor [1] - relative level of Rayleigh damping applied compared to that 
%   estimated on line 23 from diagonal of mass and stiffness matrices
%OUPUTS
%   C - same as input C but with addition of necessary damping terms on
%   diagonal
%--------------------------------------------------------------------------
if isempty(varargin)
    interface_damping_factor = 1;
else
    interface_damping_factor = varargin{1};
end

%Names of all interface elements
interface_el_names = {'ASI2D2'};

%Estimated damping level required, based on flaky calculation!
damping_level = min(sqrt(full(diag(K)) ./ full(diag(M)))) / (pi ^ 2) *  interface_damping_factor;

%Find the interface elements and associated nodes
nds_of_interface_els = els(ismember(el_typ_i, interface_el_names),:);
nds_of_interface_els = nds_of_interface_els(nds_of_interface_els > 0);
if isempty(nds_of_interface_els)
    %No action needed - there are no interface elements in model
    return
end

%Find the relevant indices for DoFs associated with all interface nodes
%in global matrix
gi = gl_lookup(nds_of_interface_els, :);
gi = gi(:);
gi = gi(gi > 0);

% gi = gl_lookup(nds_of_interface_els, 4);
% gi = gi(:);
% gi = gi(gi > 0);

%Add Rayleigh damping on the diagonal of the damping matrix for these DoFs
%proportional to mass terms
C(gi,gi) = C(gi,gi) + M(gi,gi) * damping_level;

end