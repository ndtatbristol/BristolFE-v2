clear
close all
addpath(genpath('../code'));

%Distance to faceted 2D boundary

no_pts_per_side = 11; %use <= 10 pts to see nearest points and normals
a = linspace(-10,10,no_pts_per_side);
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
if no_pts_per_side > 10
    d1 = fn_signed_dist_to_bdry([x(:), y(:)], bdry_vtcs);
else
    [d1, nearest_pts, norm_vecs] = fn_signed_dist_to_bdry([x(:), y(:)], bdry_vtcs, bdry_fcs, interior_pt);
end
toc



d1 = reshape(d1, size(x));

figure;
imagesc(a,a,d1);
hold on;
plot([bdry_vtcs(bdry_fcs(:, 1), 1), bdry_vtcs(bdry_fcs(:, 2), 1)]', ...
     [bdry_vtcs(bdry_fcs(:, 1), 2), bdry_vtcs(bdry_fcs(:, 2), 2)]', ...
    'r')
axis equal;
if numel(a) <= 10
    plot(nearest_pts(:,1), nearest_pts(:,2), 'gx');
    plot([nearest_pts(:,1), nearest_pts(:,1) + norm_vecs(:,1)]', [nearest_pts(:,2), nearest_pts(:,2) + norm_vecs(:,2)]', 'g.-');
end
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
