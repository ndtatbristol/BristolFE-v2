clear
close all
addpath(genpath('../code'));

%Some test points
a = linspace(-10,10,100);
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
d = fn_signed_dist_to_bdry([x(:), y(:), z(:)], bdry_nds, bdry_fcs, interior_pt);
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

figure;
imagesc(a,a,d(:,:,50));
axis equal;
colorbar
