clear
close all
m1 = 1;
m2 = 2;
dt = 0.01;
tpts = 10000;

c11 = 0;
c12 = -0.1;
c21 = 0.1;
c22 = 0;

inv_M = [1/m1, 0; 0, 1/m2];
C = [c11, c12; c21, c22];

A = [2 * eye(2) - dt * inv_M * C,   dt * inv_M * C - eye(2);
     eye(2),                        zeros(2)];

u = zeros(tpts,2);
u(1, 1) = 1;
for i = 3:tpts
    u(i, :) = (A(1:2,:) * [u(i-1,:)'; u(i-2,:)'])';
end

figure;plot(u)