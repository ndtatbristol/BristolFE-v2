clear
close all
addpath(genpath('../code'));

%2 Distance to faceted boundary
a = linspace(-10,10,200);
[x,y] = meshgrid(a,a);
bdry_vtcs = [
    0, 0
    1, 0
    1, 2
    0, 1] * 3;
bdry_fcs = [
    1, 2
    2, 3
    4, 3
    1, 4];
interior_pt = [2, 2];

tic
% d1 = fn_signed_dist_to_bdry([x(:), y(:)], bdry_vtcs);
% d1 = fn_signed_dist_to_bdry([x(:), y(:)], bdry_vtcs, bdry_fcs);
d1 = fn_signed_dist_to_bdry([x(:), y(:)], bdry_vtcs, bdry_fcs, interior_pt);
toc



d1 = reshape(d1, size(x));

figure;
imagesc(a,a,d1);
hold on;
plot([bdry_vtcs(bdry_fcs(:, 1), 1), bdry_vtcs(bdry_fcs(:, 2), 1)]', ...
     [bdry_vtcs(bdry_fcs(:, 1), 2), bdry_vtcs(bdry_fcs(:, 2), 2)]', ...
    'r')
axis equal;
c = caxis;
colorbar

%comparison with old function
% tic
% d1_old = fn_dist_point_to_bdry_2D([x(:), y(:)], bdry_vtcs); %old version
% toc
% d1_old = reshape(d1_old, size(x));
% figure;
% imagesc(a,a,d1_old);
% hold on;
% plot(bdry_vtcs([1:end, 1],1), bdry_vtcs([1:end, 1],2), 'r')
% title('Old version')
% axis equal;
% caxis(c);
% colorbar
