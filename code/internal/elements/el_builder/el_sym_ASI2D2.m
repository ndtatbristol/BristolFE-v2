%symbolic creation of acoustic interface 2D element
%Normal should point into fluid for consistency with abaqus: https://classes.engineering.wustl.edu/2009/spring/mase5513/abaqus/docs/v6.6/books/usb/default.htm?startat=pt10eli02.html
%which relies on node order being correct!!!!
restoredefaultpath;
addpath(genpath('../..'));

clear;
clc;

no_nds = 2; %number of nodes of element
no_dims = 2; %this is number of dimensions of element
dof_indices = [1,2,3,4]; %indices of the DOF at each node

el_type = 'ASI2D2';

%--------------------------------------------------------------------------
%Change current folder to where this file is to ensure generated element
%files are in right place (i.e. one folder up);
cd(fileparts(mfilename('fullpath')));

no_dfs = numel(dof_indices);

[loc_nd, loc_df] = meshgrid([1:no_nds], dof_indices);
loc_nd = loc_nd(:);
loc_df = loc_df(:);


%Define symbols
nds = sym('nds_%d_%d', [no_nds, no_dims], 'real');

%Calculate length, detJ, and unit normal
tmp = nds(1,:) - nds(2,:);
detJ = sqrt(sum(tmp .^ 2));
unit_norm = [tmp(2), -tmp(1)] / detJ;

K = sym(zeros(no_dfs * no_nds));

C = sym(zeros(no_dfs * no_nds));
switch el_type
    case 'ASI2D2_a'
        C([1,2,5,6],[4,8]) = -detJ / 2 * [unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2); unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2)];
        C([4,8],[1,2,5,6]) = -detJ / 2 * [unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2); unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2)].';
    case 'ASI2D2_b'
        C([1,2,5,6],[4,8]) = detJ / 2 * [unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2); unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2)];
        C([4,8],[1,2,5,6]) = detJ / 2 * [unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2); unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2)].';
    case 'ASI2D2_c'
        C([1,2,5,6],[4,8]) = -detJ / 4 * [unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2); unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2)];
        C([4,8],[1,2,5,6]) = -detJ / 4 * [unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2); unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2)].';
    otherwise %this is the right one!!
        C([1,2,5,6],[4,8]) = -detJ / 4 * [unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2); unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2)];
        C([4,8],[1,2,5,6]) = -detJ / 4 * [unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2); unit_norm(1), unit_norm(1); unit_norm(2), unit_norm(2)].';
end

M = sym(zeros(no_dfs * no_nds));

%--------------------------------------------------------------------------
fn_create_element_matrix_file(['..', filesep, 'fn_el_', el_type, '.m'], K, C, M, detJ, loc_nd, loc_df, no_dims);

fn_el_mats = str2func(['fn_el_', el_type]);

D = rand(6);
D = D + D';
rho = 1234.5;

%Test function #1 - limitied DOF, small
[el_K, el_C, el_M, loc_nd, loc_df] = fn_el_mats([0,0; cosd(30), sind(30)], [1, 2], D, rho);
disp(squeeze(el_C));

% [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_ASI2D2_a([0,0; cosd(30), sind(30)], [1, 2], D, rho);
% disp(squeeze(el_C));
% [el_K, el_C, el_M, loc_nd, loc_df] = fn_el_ASI2D2_b([0,0; cosd(30), sind(30)], [1, 2], D, rho);
% disp(squeeze(el_C));


%Test function #1 - all DOF, big
tic;
n = 100000;
[el_K2, el_C2, el_M2, loc_nd2, loc_df2] = fn_el_mats(rand(n,2), randi(n,n,3), D, rho);
toc;

