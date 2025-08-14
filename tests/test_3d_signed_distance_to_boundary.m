clear
close all
addpath(genpath('../code'));

%Some test points
a = linspace(-10,10,10);
[x,y,z] = meshgrid(a,a,a);


bdry_nds = [
    -1, 0, -1
    0, -1, -1
    1, 1, 0
    -1, -1, 2] * 3;

bdry_fcs = [
    1,2,3
    3,2,4
    1,2,4
    4,1,3];

interior_pt = [0,0,0];

tic
if numel(a) > 10;
    d = fn_signed_dist_to_bdry([x(:), y(:), z(:)], bdry_nds, bdry_fcs, interior_pt);
else
    [d, nearest_pts, norm_vecs]  = fn_signed_dist_to_bdry([x(:), y(:), z(:)], bdry_nds, bdry_fcs, interior_pt);
end
toc

d = reshape(d, size(x));
figure;
xslice = 0;
yslice = 0;
zslice = 0;
h = slice(x,y,z,d,xslice,yslice,zslice);
for i = 1 :numel(h)
    set(h(i), 'EdgeColor', 'None')
end
hold on;
patch('Faces', bdry_fcs, 'Vertices', bdry_nds,'FaceColor', 'r', 'FaceAlpha', 0.5);
axis equal
colorbar
if numel(a) <= 10
    plot3(nearest_pts(:,1), nearest_pts(:,2), nearest_pts(:,3), 'gx');
    plot3([nearest_pts(:,1), nearest_pts(:,1) + norm_vecs(:,1)]', ...
        [nearest_pts(:,2), nearest_pts(:,2) + norm_vecs(:,2)]', ...
        [nearest_pts(:,3), nearest_pts(:,3) + norm_vecs(:,3)]', ...
        'g.-');
end


% figure;
% imagesc(a,a,d(:,:,round(numel(a) / 2)));
% axis equal;
% colorbar
