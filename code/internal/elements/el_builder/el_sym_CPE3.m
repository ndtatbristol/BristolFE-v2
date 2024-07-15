%symbolic creation of mass and stiffness matrices
%for 3 node triangular linear (constant strain) elastic element CPE3
restoredefaultpath;
addpath(genpath('../..'));

clear;
clc;

%NOTE - the material stiffness matrix, D, needs to be defined in terms of 
%engineering shear strain, e.g.
%   gamma_xy = du_x/dy + du_y/dx = 2 * epsilon_xy
%not shear strain, which is defined as
%   epsilon_xy = (du_x/dy + du_y/dx) / 2 = gamma_xy / 2
%This means that D matrix will be  symmetric even for generally anisotropic
%materials.
%BE VERY CAREFUL about the definition used if transcribing stiffness matrix
%terms from textbooks etc. For example efunda [1] presents the stiffness 
%matrices in terms of epsilon_** for basic classes of material (isotropic, 
%orthotropic, transversely isotropic etc). For these classes of material in 
%the orientations presented, the last 3 columns (those related to shear 
%strains) are all zeros except on the diagonal so the matrices are symmetric
%regardless of the definition of shear strain (although there is a note - not 
%obvious - in the text below the matrix on the efunda page explaining that 
%they use epsilon_** = gamma_** / 2).

%[1] https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_isotropic.cfm

no_nds = 3; %number of nodes of element
no_dims = 2; %this is number of dimensions of element
dof_indices = [1,2,3]; %indices of the DOF at each node
no_stress = 6; %number of stress components
phys_dims = 3; %do not alter!!!

el_type = 'CPE3';
el_type = 'test';

%--------------------------------------------------------------------------
%Change current folder to where this file is to ensure generated element 
%files are in right place (i.e. one folder up);
cd(fileparts(mfilename('fullpath'))); 

%Define symbols
nds = sym('nds_%d_%d', [no_nds, no_dims], 'real'); 
rho = sym('rho');
k = sym('k');
D = sym('D_%d_%d', [no_stress, no_stress]); %stiffness matrix
q = sym('q', [1, no_dims], 'real'); %natural coordinates

%Shape functions in natural coordinates - # shape fncs = # nodes
n = sym('n', [1, no_nds], 'real'); %shape functions
n(1) = q(1);
n(2) = q(2);
n(3) = 1 - q(1) - q(2);

[detJ, invJ_times_detJ] = fn_symbolic_inv_jacobian(n, nds, q);

%Calculation of B matrix (internal strain from nodal displacements)
diff_matrix = [
    1, 0, 0
    0, 2, 0
    0, 0, 3
    0, 3, 2
    3, 0, 1
    2, 1, 0];
[B, J, N, loc_nd, loc_df] = fn_symbolic_B_matrix(diff_matrix, n, invJ_times_detJ, no_nds);

%Integrate to get K matrix
integrand = B' * D * B * J; 
z= sym('z', 'real'); %dummy axis for inegration through thickness
K = simplify(int(int(int(integrand, q(1), 0, 1 - q(2)), q(2), 0, 1), z, 0, 1));

%Make it a diagonal mass matrix - would be more logical to just share total
%mass to nodes directly rather than going through consistent mass matrix
%formualtion and then diagonalising!
M = J * simplify(int(int(N' * N, q(1), 0, 1 - q(2)), q(2), 0, 1)) * rho;
M = diag(sum(M));

%Remove rows/cols for unwanted DOF from symbolic K and M matrices
j = ismember(loc_df, dof_indices);
K = K(j, j);
M = M(j, j);
loc_df = loc_df(j);
loc_nd = loc_nd(j);

C = sym(zeros(size(K)));
%--------------------------------------------------------------------------
fn_create_element_matrix_file(['..', filesep, 'fn_el_', el_type, '.m'], K, C, M, detJ, loc_nd, loc_df, no_dims);

fn_el_mats = str2func(['fn_el_', el_type]);

D = rand(6);
D = D + D';
rho = 1234.5;

%Test function #1 - limited DOF, just one element
[el_K, el_C, el_M, loc_nd, loc_df] = fn_el_mats([0,0;1,0;0,1], [1,2,3], D, rho, [1,2]);
disp(squeeze(el_K));
disp(squeeze(el_M));
return

%Test function #1 - all DOF, big
tic;
n = 100000;
[el_K2, el_C2, el_M2, loc_nd2, loc_df2] = fn_el_mats(rand(n,2), randi(n,n,3), D, rho);
toc;

