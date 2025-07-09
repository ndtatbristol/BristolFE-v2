clear
% close all
addpath(genpath('../code'));

%Some test points
a = linspace(-10,10,100);
[x,y,z] = meshgrid(a,a,a);

if 0
    %1. Test distance to plane function for multiple pts
    a = [3,0,0];
    b = [0,3,0];
    c = [5,5,5];

    [d, alpha, beta] = fn_dist_point_to_plane([x(:), y(:), z(:)], a + c, b + c, c);
    d = reshape(d, size(x));
    alpha = reshape(alpha, size(x));
    beta = reshape(beta, size(x));

    tmp = double((alpha >= 0) & (beta >= 0) & ((alpha + beta) <= 1));

    figure;
    xslice = [1, 5, 9];
    yslice = 5;
    zslice = 5;
    slice(x,y,z,tmp,xslice,yslice,zslice)
end

if 1
    %2 Distance to faceted boundary
    bdry_nds = [
        -5, 1, 1
        1, 5, 1
        1,-5, -2
        1, 1, 5] / 5 * 8;
    bdry_nds = [
        -5, 0, 0
        0, 5, 0
        0,-5, 0
        0, 0, 5] / 5 * 8;
    bdry_fcs = [1,2,3
        2,3,4];
    tic
    d = fn_dist_point_to_bdry_3D([x(:), y(:), z(:)], bdry_nds, bdry_fcs);
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
end