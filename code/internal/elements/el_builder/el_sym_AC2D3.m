%symbolic creation of mass and stiffness matrices
%for 3 node triangular acoustic element AC2D3
restoredefaultpath;
addpath(genpath('../FE functions'));

clear;
clc;

no_nds = 3; %number of nodes of element
no_dims = 2; %this is number of dimensions of element
dof_indices = 1; %indices of the DOF at each node
phys_dims = 3; %do not alter!!!

pressure_dof_i = 4;

el_type = 'AC2D3_v2';
%--------------------------------------------------------------------------

%Define symbols
nds = sym('nds_%d_%d', [no_nds, no_dims], 'real'); 
rho = sym('rho');
k = sym('k');
D = sym('D');
q = sym('q', [1, no_dims], 'real'); %natural coordinates

%Shape functions in natural coordinates - # shape fncs = # nodes
n = sym('n', [1, no_nds], 'real'); %shape functions
n(1) = q(1);
n(2) = q(2);
n(3) = 1 - q(1) - q(2);

[detJ, invJ_times_detJ] = fn_symbolic_inv_jacobian(n, nds, q);

%SHould be correct up to here - same as CST elastic element

%Calculation of B matrix (internal strain from nodal displacements)
diff_matrix = [
    1
    2
    3];
[B, J, N, loc_nd, loc_df] = fn_symbolic_B_matrix(diff_matrix, n, invJ_times_detJ, no_nds);

%Integrate to get K matrix
integrand = B' * B * J; 
z= sym('z', 'real'); %dummy axis for inegration through thickness
K = simplify(int(int(int(integrand, q(1), 0, 1 - q(2)), q(2), 0, 1), z, 0, 1));

%Make it a diagonal mass matrix - would be more logical to just share total
%mass to nodes directly rather than going through consistent mass matrix
%formualtion and then diagonalising!
M = J * simplify(int(int(N' * N, q(1), 0, 1 - q(2)), q(2), 0, 1));
M = diag(sum(M));

%The critical lines that determine what the field variable actually is!
K = -K / rho;
M = -M / D;

%Remove rows/cols for unwanted DOF from symbolic K and M matrices
j = ismember(loc_df, dof_indices);
K = K(j, j);
M = M(j, j);
loc_df = loc_df(j) * pressure_dof_i;
loc_nd = loc_nd(j);

C = sym(zeros(size(K)));
%--------------------------------------------------------------------------
%Minus signs before K, M and C addded 29/3/24 by PDW to avoid instablity in
%fluid solid interfaces (this is what enables use of symmetric coupling
%matrix for interface elements).
fn_create_element_matrix_file(['fn_el_', el_type, '.m'], -K, -C, -M, detJ, loc_nd, loc_df, no_dims);%, '%%Velocity squared', 'c2 = D / rho;');

fn_el_mats = str2func(['fn_el_', el_type]);

D = rand(1);
rho = 1234.5;

%Test function #1 - limitied DOF, small
[el_K, el_M, loc_nd, loc_df] = fn_el_mats(rand(10,2), randi(10,5,3), D, rho, [1,2]);

%Test function #1 - all DOF, big
tic;
n = 100000;
[el_K2, el_M2, loc_nd2, loc_df2] = fn_el_mats(rand(n,2), randi(n,n,3), D, rho);
toc;

fname = ['fn_el_', el_type, '.m'];
copyfile(fname, '../FE functions');
fprintf([fname, '\n']);