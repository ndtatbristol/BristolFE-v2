function max_EV = fn_max_EV_of_simple_interface(k1, m1, c1, k2, m2, c2, dt, coupling_eff)

ndf = 2;

K = [k1,    0
     0,     k2];

M = [m1,    0
     0,     m2];

inv_M = [1 / m1,    0
         0,         1 / m2];

C = [c1,    coupling_eff
     coupling_eff,     c2];

%Prop matrix based on last two disps
B1 = 2 * eye(ndf) - dt * (inv_M * C) - dt ^ 2 * (inv_M * K);
B2 = dt * (inv_M * C) - eye(ndf);
A = [B1,        B2
     eye(ndf),  zeros(ndf)];

%Prop matrix based on last disp and last vel
% A = [eye(ndf) - dt ^ 2 * inv_M * K,     dt * (eye(ndf) - dt * inv_M * C)
%               - dt     * inv_M * K,           eye(ndf) - dt * inv_M * C];

%Velocity at the current step
% B1 = 2 * speye(ndf) - dt ^ 2 * inv_M * K;
% B2 =    -speye(ndf) + dt / 2 * inv_M * C;
% B3 =     speye(ndf) + dt / 2 * inv_M * C;
% A = [B3 \ B1,   B3 \ B2
%     eye(ndf), zeros(ndf)];

% max_EV = max(abs(eig(A)));
max_EV = abs(eigs(A, 1));

end